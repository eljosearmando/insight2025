% Glucose–insulin simulation with subcutaneous insulin absorption
% and a simple meal glucose appearance model.
clear; clc;

%% --------------------- USER SETTINGS ---------------------
T_end_min     = 6*60;        % total sim time [min]
dt_plot       = 0.1;         % plot resolution [min]

% Meal (carbohydrates)
meal_time_min = 30;          % meal time [min]
meal_size_g   = 60;          % grams of carbs

% Insulin dosing (rapid-acting SC bolus + optional basal)
bolus_time_min  = 20;        % bolus time [min]
bolus_units     = 6;         % units
basal_U_per_hr  = 0.8;       % units/hour continuous basal (set 0 for none)

% Initial conditions (fasting steady-ish)
Gb   = 100;    % basal glucose [mg/dL]
Ib   = 10;     % basal plasma insulin [mU/L]
X0   = 0;      % insulin action compartment

%% ------------------- MODEL PARAMETERS --------------------
Vg = 12;                   % glucose effective volume [dL]
VI = 12;                   % insulin plasma distribution volume [L]
n  = 0.14;                 % insulin clearance rate [min^-1]
p1 = 0.028;                % glucose effectiveness [min^-1]
p2 = 0.025;                % insulin action decay [min^-1]
p3 = 1e-4;                 % insulin action sensitivity [(min^-2)/(mU/L)]

ka1 = 0.035;               % SC depot -> SC transit [min^-1]
ka2 = 0.03;                % SC transit -> plasma [min^-1]
U_to_mU = 6000;            % 1 U = 6000 mU

% Meal appearance (two-compartment gut model)
f_bio   = 0.9;             % fraction of carbs that reach plasma
tau_gut = 25;              % gut time constant [min]
k_g1    = 1/tau_gut;       % gut emptying rate [min^-1]
k_g2    = 1/(tau_gut*1.3); % slower transit
g_to_mg = 1000;            % convert g carbs -> mg glucose

% Basal insulin infusion (convert U/hr to mU/min into SC depot)
basal_mU_per_min = basal_U_per_hr * U_to_mU / 60;

%% --------------------- SIMULATION ------------------------
% State vector y = [G, X, I, S1, S2, Q1, Q2]'
y0 = [Gb; X0; Ib; 0; 0; 0; 0];
tspan = [0 T_end_min];

params = struct('Gb',Gb,'Ib',Ib,'Vg',Vg,'VI',VI,'n',n,'p1',p1,'p2',p2,'p3',p3, ...
   'ka1',ka1,'ka2',ka2,'U_to_mU',U_to_mU,'g_to_mg',g_to_mg, ...
   'f_bio',f_bio,'k_g1',k_g1,'k_g2',k_g2, ...
   'basal_mU_per_min',basal_mU_per_min, ...
   'bolus_time',bolus_time_min,'bolus_units',bolus_units, ...
   'meal_time',meal_time_min,'meal_size_g',meal_size_g, ...
   ... %%% CHANGED: smooth meal absorption parameters
   'meal_mg', meal_size_g * g_to_mg * f_bio, ...   % total mg to appear
   'tau_ingest', 25, ...                            % min, controls ingestion spread
   ... %%% CHANGED: exercise parameters (time-windowed effects)
   'ex_start',  45, ...       % start of exercise [min]
   'ex_dur',    45, ...       % duration [min]
   'ex_int',   0.7, ...       % intensity 0..1 (moderate)
   'ex_boost_p1', 0.50, ...   % +50% p1 during exercise (scaled by intensity)
   'ex_boost_p3', 0.50, ...   % +50% p3 during exercise (scaled by intensity)
   'k_ex', 0.01, ...          % optional direct sink (1/min) * intensity
   'post_gain_p3', 0.30, ...  % +30% p3 immediately after exercise
   'post_tau_min', 120 ...    % decay time constant [min] for after-effect
   );

opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) gi_ode(t,y,params), tspan, y0, opts);

% Interpolate for smooth plots
tq = (0:dt_plot:T_end_min).';
yq = interp1(t,y,tq,'pchip');

G  = yq(:,1); X = yq(:,2); I = yq(:,3);
S1 = yq(:,4); S2 = yq(:,5); Q1 = yq(:,6); Q2 = yq(:,7);

% Compute inputs/appearance rates for plotting (vectorized)
[u_sc, u_plasma_mU_per_min, Ra_mg_per_min] = gi_inputs(tq, yq, params);

%% ----------------------- PLOTTING ------------------------
figure('Color','w');

subplot(3,1,1)
plot(tq,G,'LineWidth',1.8); grid on
ylabel('Glucose (mg/dL)'); xlabel('Time (min)')
title('Plasma Glucose')
xline(meal_time_min,'--','Meal'); xline(bolus_time_min,'--','Bolus');
xline(params.ex_start,'k:','Ex start'); xline(params.ex_start+params.ex_dur,'k:','Ex end'); %%% CHANGED (optional markers)

subplot(3,1,2)
plot(tq,I,'LineWidth',1.8); grid on
ylabel('Insulin (mU/L)'); xlabel('Time (min)')
title('Plasma Insulin')
xline(meal_time_min,'--'); xline(bolus_time_min,'--');
xline(params.ex_start,'k:'); xline(params.ex_start+params.ex_dur,'k:'); %%% CHANGED

subplot(3,1,3)
plot(tq, Ra_mg_per_min/params.Vg, 'LineWidth',1.8); hold on; grid on
plot(tq, u_plasma_mU_per_min/params.VI, 'LineWidth',1.8);
ylabel('Rates (per dL or per L)'); xlabel('Time (min)')
legend({'Glucose appearance / V_g','Insulin appearance / V_I'},'Location','best')
title('Appearance Rates (normalized by distribution volumes)')
xline(meal_time_min,'--'); xline(bolus_time_min,'--');
xline(params.ex_start,'k:'); xline(params.ex_start+params.ex_dur,'k:'); %%% CHANGED

sgtitle('Insulin Activity Simulation (Meal + SC Bolus + Basal + Exercise)') %%% CHANGED

%% --------------------- METRICS DISPLAY -------------------
G_min = min(G); G_max = max(G); G_end = G(end);
fprintf('--- Key Glucose Metrics ---\n');
fprintf('Min glucose: %.1f mg/dL\n', G_min);
fprintf('Max glucose: %.1f mg/dL\n', G_max);
fprintf('End glucose: %.1f mg/dL at t = %.0f min\n', G_end, T_end_min);

%% ==================== LOCAL FUNCTIONS ====================
function dy = gi_ode(t, y, P)
   % ----- Unpack scalar state vector (ode45 passes 7x1) -----
   G  = y(1);  % mg/dL
   X  = y(2);  % 1/min
   I  = y(3);  % mU/L
   S1 = y(4);  % mU
   S2 = y(5);  % mU
   Q1 = y(6);  % mg
   Q2 = y(7);  % mg

   % ---------- Inputs ----------
   % Basal SC insulin (mU/min)
   u_basal = P.basal_mU_per_min;

   % Bolus insulin as a narrow pulse (mU/min)
   bolus_area_mU = P.bolus_units * P.U_to_mU;
   pulse_w = 0.01; % min
   if abs(t - P.bolus_time) <= pulse_w/2
       u_bolus = bolus_area_mU / pulse_w;
   else
       u_bolus = 0;
   end
   u_sc = u_basal + u_bolus; % mU/min into S1

   %%%% CHANGED: Meal ingestion as smooth exponential into Q1 (no impulse)
   if t >= P.meal_time
       u_meal = (P.meal_mg / P.tau_ingest) * exp(-(t - P.meal_time)/P.tau_ingest); % mg/min
   else
       u_meal = 0;
   end

   % ---------- SC insulin dynamics ----------
   dS1dt = -P.ka1*S1 + u_sc;          % mU/min
   dS2dt =  P.ka1*S1 - P.ka2*S2;      % mU/min
   u_plasma_mU_per_min = P.ka2*S2;    % mU/min into plasma

   % ---------- Exercise effects (time-windowed)  %%% CHANGED
   in_ex   = (t >= P.ex_start) && (t <= P.ex_start + P.ex_dur);
   post_ex = (t >  P.ex_start + P.ex_dur);
   dt_post = t - (P.ex_start + P.ex_dur);

   mult_p1 = 1.0;   mult_p3 = 1.0;   sink_ex = 0;
   if in_ex
       mult_p1 = mult_p1 * (1 + P.ex_boost_p1 * P.ex_int);
       mult_p3 = mult_p3 * (1 + P.ex_boost_p3 * P.ex_int);
       sink_ex = P.k_ex * P.ex_int;   % extra disposal ~ exercise
   elseif post_ex
       mult_p3 = mult_p3 * (1 + P.post_gain_p3 * exp(-dt_post / P.post_tau_min));
   end
   p1_eff = P.p1 * mult_p1;
   p3_eff = P.p3 * mult_p3;

   % ---------- Plasma insulin kinetics ----------
   dIdt = -P.n*(I - P.Ib) + (u_plasma_mU_per_min / P.VI); % mU/L/min

   % ---------- Gut glucose appearance ----------
   dQ1dt = -P.k_g1*Q1 + u_meal;     % mg/min
   dQ2dt =  P.k_g1*Q1 - P.k_g2*Q2;  % mg/min
   Ra_mg_per_min = P.k_g2*Q2;       % mg/min
   Ra_over_Vg = Ra_mg_per_min / P.Vg; % mg/dL/min

   % ---------- Insulin action ----------
   dXdt = -P.p2*X + p3_eff*(I - P.Ib); % 1/min  %%% CHANGED (uses p3_eff)

   % ---------- Glucose dynamics ----------
   dGdt = -(X + p1_eff)*(G - P.Gb) + Ra_over_Vg - sink_ex*(G - P.Gb); %%% CHANGED

   % Pack derivatives (column vector 7x1)
   dy = [dGdt; dXdt; dIdt; dS1dt; dS2dt; dQ1dt; dQ2dt];
end

function [u_sc, u_plasma_mU_per_min, Ra_mg_per_min] = gi_inputs(t, Y, P)
   % Vectorized helper for plotting inputs/appearance over many times t.
   % Y is N x 7 matrix of states at those times.
   S2 = Y(:,5);
   Q2 = Y(:,7);

   % Basal SC insulin (vector)
   u_basal = P.basal_mU_per_min * ones(size(t));

   % Bolus pulse over grid (vector)
   bolus_area_mU = P.bolus_units * P.U_to_mU;
   pulse_w = 0.01; % min (match ODE)
   u_bolus = zeros(size(t));
   u_bolus(abs(t - P.bolus_time) <= pulse_w/2) = bolus_area_mU / pulse_w;

   u_sc = u_basal + u_bolus;                 % mU/min into S1
   u_plasma_mU_per_min = P.ka2 .* S2;        % mU/min into plasma
   Ra_mg_per_min = P.k_g2 .* Q2;             % mg/min into plasma
end
