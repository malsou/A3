function y = eval_engine(des, op, opts)
%EVAL_ENGINE Low-order coupled cycle + cooling evaluation for DOE screening.
% Required design fields (des):
%   BPR, FPR, LPR, HPR, Tt4, nozzle, cooling, porosity, tTBC, kSt
% Required operating-point fields (op struct array):
%   alt_ft, M0, Freq
% Optional opts fields are provided by engine_defaults().
%
% Output fields retained for A3 pipeline compatibility:
%   y.TSFC_cruise, y.phi_crit, y.MT_worst, y.MF_worst, y.J, y.P

if nargin < 3 || isempty(opts)
    opts = engine_defaults();
else
    opts = engine_defaults(opts);
end

OPR = des.LPR * des.HPR;

% Coupled compressor efficiency from map-position proxy
eta_pc = opts.eta_pc_peak ...
    - opts.eta_pc_a*(log(OPR/opts.OPR0))^2 ...
    - opts.eta_pc_b*(opts.SMc-opts.SMc0)^2;
eta_pc = clamp(eta_pc, opts.eta_pc_bounds(1), opts.eta_pc_bounds(2));

isMixed = strcmpi(des.nozzle, "mixed");
isEffusion = strcmpi(des.cooling, "effusion");
nozFactor_TSFC = 1.0 + opts.nozFactorMixed*isMixed;
coolFactor_eff = opts.coolEffEffusion*isEffusion;

ST = opts.ST_base ...
    + opts.ST_Tt4_gain*(des.Tt4 - 1500) ...
    + opts.ST_BPR_gain*(des.BPR - 8) ...
    + opts.ST_OPR_gain*(OPR - 35);
ST = max(ST, opts.ST_min);

atm0 = atmos_isa(op(1).alt_ft);
lapseSLS = thrust_lapse(atm0, op(1).M0);
mdot_core = op(1).Freq / (ST * lapseSLS);
mdot_core = max(mdot_core, opts.mdot_core_min);

eta_th = 0.30 + 0.08*log(OPR)/log(40) + 0.06*((eta_pc - 0.89)/0.03);
eta_th = clamp(eta_th, opts.eta_th_bounds(1), opts.eta_th_bounds(2));

eta_p = 0.55 + 0.22*((des.BPR - 4)/8) - 0.10*((des.FPR - 1.60)/0.30);
eta_p = clamp(eta_p, opts.eta_p_bounds(1), opts.eta_p_bounds(2));

n = numel(op);
TSFC = zeros(n,1);
phi  = zeros(n,1);
Tmax = zeros(n,1);
Fnet = zeros(n,1);
Freq = zeros(n,1);

for i = 1:n
    atm = atmos_isa(op(i).alt_ft);
    lapse = thrust_lapse(atm, op(i).M0);

    [phi(i), Tmax(i)] = cooling_requirement(des, OPR, atm, op(i).M0, coolFactor_eff, opts);

    % Coupled turbine efficiency penalty from cooling requirement
    eta_pt = opts.eta_pt_peak - opts.eta_pt_penalty_c*phi(i);
    eta_pt = clamp(eta_pt, opts.eta_pt_bounds(1), opts.eta_pt_bounds(2));

    bleedFactor = max(1 - 2.0*phi(i), 0.7);
    mixPenalty = 1 + 6.0*phi(i);

    Fnet(i) = mdot_core * ST * lapse * bleedFactor;
    TSFC(i) = (1.0/(eta_th*eta_p)) * (0.93/eta_pt) ...
        * (1 + 0.10*((des.Tt4 - 1500)/300)) ...
        * nozFactor_TSFC * mixPenalty;

    Freq(i) = op(i).Freq;
end

MT = opts.T_allow - Tmax;
MF = (Fnet - Freq)./Freq;

y.MT_worst = min(MT);
y.MF_worst = min(MF);
y.phi_crit = max(phi);

[~, iCruise] = max([op.alt_ft]);
y.TSFC_cruise = TSFC(iCruise);

if n >= 2
    [~, idxSort] = sort([op.alt_ft]);
    iClimb = idxSort(max(1, n-1));
else
    iClimb = iCruise;
end

y.J = opts.J_weights(1)*y.TSFC_cruise + opts.J_weights(2)*TSFC(iClimb);
y.P = max(0,-y.MF_worst) + opts.beta*max(0,-y.MT_worst/opts.Tscale) + opts.gamma*opts.P_solver;
end

function lapse = thrust_lapse(atm, M0)
% Density and Mach-based thrust lapse proxy.
rho0 = 1.225;
sigma = atm.rho / rho0;
lapse = sigma * (1 - 0.15*M0);
lapse = max(lapse, 0.2);
end

function [phi, Tmax] = cooling_requirement(des, OPR, atm, M0, eff_bonus, opts)
% Cooling requirement model using heat-transfer-strength scaling and
% film/effusion effectiveness concept motivated by turbine cooling notes.

porN = (des.porosity - 0.002)/(0.020-0.002);
CdN  = (opts.Cd - 0.60)/(0.95-0.60);
angN = sin(deg2rad(opts.alpha));
tN   = (des.tTBC - 100)/(400-100);

% Effectiveness proxy bounded to realistic range
eff = 0.35 + 0.25*porN + 0.15*CdN + 0.10*angN + 0.15*tN + eff_bonus;
eff = clamp(eff * des.kSt, 0.20, 0.90);

% External gas temperature proxy with altitude/Mach effects
Tg = des.Tt4 - 120*(atm.h_m/11000) + 20*M0;

% Stanton / HTC scaling: kSt multiplies transfer intensity
hScale = des.kSt * (1 + 0.25*((OPR - 20)/35) + 0.05*M0);
hScale = max(hScale, 0.4);

% Required coolant fraction to enforce Tmax <= T_allow
phi_req = (Tg - opts.T_allow - 0.35*(des.tTBC - 100)) ...
    / (900*(0.6 + 0.4*eff)*hScale);
phi = clamp(phi_req + 0.005, opts.phi_bounds(1), opts.phi_bounds(2));

Tmax = 980 + 0.55*(des.Tt4 - 1450) ...
    - 900*phi*(0.6 + 0.4*eff)*hScale ...
    - 0.35*(des.tTBC - 100);
Tmax = clamp(Tmax, opts.Tmetal_floor, opts.Tmetal_ceil);
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end
