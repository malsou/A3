function y = evaluate_engine_model_proxy(x)
% Proxy coupled cycle + cooling model.
% Inputs x: struct with fields:
% BPR,FPR,LPR,HPR,Tt4,SMc,kSt,porosity,Cd,alpha,tTBC,nozzle,cooling

% ---- Operating points (edit as desired) ----
pts = operating_points();

% ---- Derived quantities ----
OPR = x.LPR * x.HPR;

% ---- Coupled component efficiencies ----
% Rationale: compressor/turbine efficiency are not independent DOE knobs.
% They vary with cycle pressure ratio and cooling demand to preserve coupling.
% Compressor map-coupling constants (chosen to keep eta_pc in realistic range):
% eta_pc = eta_pc_peak - a*(log(OPR/OPR0))^2 - b*(SMc-SMc0)^2, clamped [0.80, 0.92]
eta_pc_peak = 0.91; OPR0 = 35; SMc0 = 0.15; a = 0.090; b = 1.80;
eta_pc = eta_pc_peak - a*(log(OPR/OPR0))^2 - b*(x.SMc-SMc0)^2;
eta_pc = min(max(eta_pc, 0.80), 0.92);

% ---- Categorical multipliers ----
isMixed = strcmpi(x.nozzle, "mixed");
isEffusion = strcmpi(x.cooling, "effusion");

nozFactor_TSFC = 1.00 - 0.015*isMixed;     % mixed slightly better
coolFactor_eff = 0.08*isEffusion;          % effusion slightly higher eff in this proxy

% ---- Cycle proxy (baseline sizing) ----
% Specific thrust proxy (N per (kg/s) of core flow)
ST = 240 ...
    + 120*((x.Tt4 - 1500)/300) ...
    - 28*((x.BPR - 8)/4) ...
    - 18*((OPR - 35)/15);
ST = max(ST, 80);

% Choose core mass flow to meet takeoff requirement at SLS (proxy sizing)
lapseSLS = thrust_lapse(pts(1).h_m, pts(1).M);
mdot_core = pts(1).F_req_N / (ST * lapseSLS);
mdot_core = max(mdot_core, 50); % avoid absurdly small flow

% Efficiency proxies (bounded)
eta_th = 0.30 ...
    + 0.08*log(OPR)/log(40) ...
    + 0.06*((eta_pc - 0.89)/0.03);
eta_th = min(max(eta_th, 0.20), 0.55);

eta_p = 0.55 ...
    + 0.22*((x.BPR - 4)/8) ...
    - 0.10*((x.FPR - 1.60)/0.30);
eta_p = min(max(eta_p, 0.40), 0.80);

% ---- Evaluate each operating point ----
TSFC = zeros(numel(pts),1);
phi  = zeros(numel(pts),1);
Tmax = zeros(numel(pts),1);
Fnet = zeros(numel(pts),1);

for i = 1:numel(pts)
    lapse = thrust_lapse(pts(i).h_m, pts(i).M);

    % Cooling proxy -> phi and Tmax (kSt scales heat-transfer effectiveness)
    [phi(i), Tmax(i)] = cooling_proxy(x, OPR, pts(i), coolFactor_eff);

    % Turbine efficiency is coupled to cooling requirement (higher phi -> lower eta_pt).
    c = 1.4; eta_pt_peak = 0.94;
    eta_pt = eta_pt_peak - c*phi(i);
    eta_pt = min(max(eta_pt, 0.82), 0.94);

    % Cooling penalties: bleed reduces thrust, mixing increases TSFC
    bleedFactor = 1 - 2.0*phi(i);                % simple bleed penalty
    bleedFactor = max(bleedFactor, 0.7);

    mixPenalty = 1 + 6.0*phi(i);                 % TSFC penalty
    % Net thrust
    Fnet(i) = mdot_core * ST * lapse * bleedFactor;

    % TSFC proxy (relative scale; consistent across samples)
    TSFC(i) = (1.0/(eta_th*eta_p)) * (0.93/eta_pt) ...
        * (1 + 0.10*((x.Tt4 - 1500)/300)) ...
        * nozFactor_TSFC ...
        * mixPenalty;
end

% ---- Responses / margins ----
T_allow = 1100; % K (edit)
MT = T_allow - Tmax;                 % margin (K)
MF = (Fnet - [pts.F_req_N]') ./ [pts.F_req_N]';  % margin (-)

% Worst-case across points
y.MT_worst = min(MT);
y.MF_worst = min(MF);

% Define "critical" coolant fraction as max across points
y.phi_crit = max(phi);

% Cruise TSFC (assume last point is cruise)
y.TSFC_cruise = TSFC(end);

% Optional objective J (weighted)
w_cruise = 0.75; w_climb = 0.25;
y.J = w_cruise*TSFC(end) + w_climb*TSFC(2);

% Penalty metric (continuous)
Tscale = 100; % K
P_solver = 0; % set to >0 if your real solver fails
beta = 1.0; gamma = 10.0;
y.P = max(0,-y.MF_worst) + beta*max(0, -y.MT_worst/Tscale) + gamma*P_solver;
end

function pts = operating_points()
% Simple 3-pt mission anchors (edit numbers later)
pts(1).name = "SLS_takeoff"; pts(1).h_m = 0;     pts(1).M = 0.0;  pts(1).F_req_N = 120e3;
pts(2).name = "climb";       pts(2).h_m = 6000;  pts(2).M = 0.4;  pts(2).F_req_N = 60e3;
pts(3).name = "cruise";      pts(3).h_m = 11000; pts(3).M = 0.78; pts(3).F_req_N = 25e3;
end

function lapse = thrust_lapse(h_m, M)
% Very rough thrust lapse (placeholder)
sigma = exp(-h_m/7000);      % density ratio proxy
lapse = sigma*(1 - 0.15*M);  % mild Mach effect
lapse = max(lapse, 0.2);
end

function [phi, Tmax] = cooling_proxy(x, OPR, pt, eff_bonus)
% Correlation-ish proxy:
% phi increases with Tt4 and OPR; decreases with porosity/Cd/tTBC/angle.
% kSt multiplies convective transfer/effectiveness (higher kSt -> lower phi/Tmax).
Tt4 = x.Tt4;

por = x.porosity;                       % 0.002-0.020
porN = (por - 0.002)/(0.020-0.002);     % 0-1
CdN  = (x.Cd - 0.60)/(0.95-0.60);       % 0-1
angN = sin(deg2rad(x.alpha));           % 0-1-ish
tN   = (x.tTBC - 100)/(400-100);        % 0-1

% Effectiveness proxy (0.2 to 0.9)
eff = 0.35 + 0.25*porN + 0.15*CdN + 0.10*angN + 0.15*tN + eff_bonus;
eff = eff .* x.kSt;
eff = min(max(eff, 0.20), 0.90);

% Thermal load proxy
load = (Tt4 - 1350)/550;         % ~0-1
load = min(max(load, 0), 2);
oprFactor = 1 + 0.25*((OPR - 20)/35) + 0.05*(pt.M);
oprFactor = max(oprFactor, 0.7);

% Coolant fraction (0 to ~0.08)
phi = 0.005 + 0.03*load*oprFactor/(0.4 + 0.6*eff);
phi = min(max(phi, 0.0), 0.08);

% Metal temperature proxy (K): increases with Tt4, decreases with cooling and TBC
Tmax = 980 + 0.55*(Tt4 - 1450) - 900*phi*(0.6 + 0.4*eff) - 0.35*(x.tTBC - 100);
Tmax = min(max(Tmax, 750), 1400);
end
