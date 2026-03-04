function y = eval_engine(des, op, opts)
%EVAL_ENGINE Coupled cycle + cooling evaluator for AEROSP 568 DOE screening.
% des required: BPR,FPR,LPR,HPR,Tt4,nozzle,cooling,porosity,tTBC,kSt
% op(i) required: alt_ft,M0,Freq
% opts optional: use engine_defaults() values + overrides

if nargin < 3 || isempty(opts)
    opts = engine_defaults();
else
    opts = engine_defaults(opts);
end

n = numel(op);
OPR = des.LPR * des.HPR;

% Coupled compressor efficiency (not an independent DOE variable)
eta_pc = opts.eta_comp_peak ...
    - opts.eta_pc_a*(log(OPR/opts.OPR0))^2 ...
    - opts.eta_pc_b*(opts.SMc - opts.SMc0)^2;
eta_pc = clamp(eta_pc, opts.eta_pc_bounds(1), opts.eta_pc_bounds(2));

% Categorical switches
isMixed = strcmpi(des.nozzle, "mixed");
isEffusion = strcmpi(des.cooling, "effusion");
coolEffBonus = opts.coolEffEffusion * isEffusion;
nozFactor_TSFC = 1.0 + opts.nozFactorMixed*isMixed;

phi = zeros(n,1);
Tmax = zeros(n,1);
Fspec = zeros(n,1);
TSFC = zeros(n,1);
Fnet = zeros(n,1);
failMetric = zeros(n,1);

for i = 1:n
    atm = atmos_isa(op(i).alt_ft);
    [phi(i), Tmax(i)] = cooling_requirement(des, OPR, atm, op(i).M0, coolEffBonus, opts);

    % Coupled turbine efficiency from coolant fraction
    eta_pt = opts.eta_t_uncooled - opts.eta_pt_penalty_c * phi(i);
    eta_pt = clamp(eta_pt, opts.eta_pt_bounds(1), opts.eta_pt_bounds(2));

    cyc = run_cycle_point(des, atm, op(i).M0, eta_pc, eta_pt, phi(i), opts);

    Fspec(i) = cyc.Fspec_net; % N per (kg/s core)
    failMetric(i) = cyc.fail_metric;
    TSFC_si = cyc.TSFC * nozFactor_TSFC;
    TSFC(i) = convert_tsfc_units(TSFC_si, opts.TSFC_unit);
end

% Engine sizing by SLS required thrust (keep same core across mission points)
mdot_core = op(1).Freq / soft_floor(Fspec(1), 1e-3, 20);
mdot_core = soft_floor(mdot_core, 1.0, 20);

for i = 1:n
    Fnet(i) = mdot_core * Fspec(i);
end

MF = (Fnet - [op.Freq]')./[op.Freq]';
MT = opts.T_allow - Tmax;

y.TSFC_cruise = TSFC(end);
y.phi_crit = max(phi);
y.MT_worst = min(MT);
y.MF_worst = min(MF);

y.J = opts.J_weights(1)*TSFC(end) + opts.J_weights(2)*TSFC(max(1,n-1));
y.P = max(0,-y.MF_worst) + opts.beta*max(0,-y.MT_worst/opts.Tscale) ...
    + opts.gamma*(opts.P_solver + mean(failMetric));
end

function cyc = run_cycle_point(des, atm, M0, eta_pc, eta_pt, phi, opts)
% 1-D turbofan cycle model inspired by AEROSP530 HW6 flow station logic
% (see repo refs: eq_sheet.pdf for cycle relations).

gc = opts.gamma_cold; gh = opts.gamma_hot;
cpc = opts.cp_cold; cph = opts.cp_hot; R = opts.R;
fail_metric = 0;

Ts0 = atm.T; Ps0 = atm.P;
Tt0 = Ts0*(1 + (gc-1)/2*M0^2);
Pt0 = Ps0*(1 + (gc-1)/2*M0^2)^(gc/(gc-1));
a0  = sqrt(gc*R*Ts0);
V0  = M0*a0;

% Inlet
Tt1 = Tt0;
Pt1 = opts.eta_r * Pt0;

% Fan: 1 -> 13
Pt13 = des.FPR * Pt1;
Tt13i = Tt1*(Pt13/Pt1)^((gc-1)/gc);
Tt13 = Tt1 + (Tt13i - Tt1)/opts.eta_fan;

% Compressor: 2 -> 3 (core)
Pt2 = Pt13; Tt2 = Tt13;
Pt3 = OPR_from_des(des) * Pt2;
Tt3i = Tt2*(Pt3/Pt2)^((gc-1)/gc);
Tt3 = Tt2 + (Tt3i - Tt2)/eta_pc;

% Work per core flow
Wfan = (1 + des.BPR)*cpc*(Tt13 - Tt1);
Wcomp = cpc*(Tt3 - Tt2);

% Bleeds / cooling split
fr_cool = clamp(opts.fr_cool_base + phi, 0.0, 0.15);
m31_raw = 1 - (opts.fr_bleed + opts.fr_leak + fr_cool);
if m31_raw <= 0
    fail_metric = fail_metric + abs(m31_raw);
end
m31 = soft_floor(m31_raw, 1e-3, 20);

% Combustor 3.1 -> 4
Pt4 = Pt3*(1 - opts.dPt_comb);
den_f = opts.eta_b*opts.FHV - cph*des.Tt4;
if den_f <= 1e3
    fail_metric = fail_metric + (1e3 - den_f)/1e3;
end
f = cph*(des.Tt4 - Tt3) / soft_floor(den_f, 1e3, 10);
m4 = m31*(1+f);

% HPT drives compressor
Tt49a = des.Tt4 - Wcomp/soft_floor(m4*cph, 1e-3, 20);
Tt49i = des.Tt4 - (des.Tt4 - Tt49a)/eta_pt;
Pt49 = Pt4*(Tt49i/des.Tt4)^(gh/(gh-1));

% Cooling mix 4.9 + cool -> 4.95
mcool = fr_cool;
Tt_cool = Tt3;
m495 = m4 + mcool;
Tt495 = (m4*cph*Tt49a + mcool*cpc*Tt_cool)/soft_floor(m495*cph, 1e-3, 20);
Pt495 = Pt49;

% LPT drives fan
Tt50a = Tt495 - Wfan/soft_floor(m495*cph, 1e-3, 20);
Tt50i = Tt495 - (Tt495 - Tt50a)/eta_pt;
Pt50 = Pt495*(Tt50i/Tt495)^(gh/(gh-1));

% Nozzles with choke-check + fully expanded exit (pe=Pa)
V19 = nozzle_velocity(Pt13, Tt13, Ps0, gc, cpc, opts.Cv19);
V9  = nozzle_velocity(Pt50, Tt50a, Ps0, gh, cph, opts.Cv9);

m13 = des.BPR;
m9 = m495;

Fspec = (m13*V19 + m9*V9) - (1+des.BPR)*V0;  % N/(kg/s)
if Fspec < 0
    fail_metric = fail_metric + (-Fspec)/(1+abs(Fspec));
end
Fspec = soft_floor(Fspec, 1e-3, 20);

mf_over_m2 = f*m31;
TSFC = mf_over_m2 / soft_floor(Fspec, 1e-3, 20); % kg/(N-s)

cyc = struct('Fspec_net',Fspec,'TSFC',TSFC,'Tt3',Tt3,'Tt495',Tt495,'fail_metric',fail_metric);
end

function V = nozzle_velocity(Pt, Tt, Pa, gamma, cp, Cv)
% Throat choking check, then fully expanded exit pe=Pa.
PRcrit = ((gamma+1)/2)^(gamma/(gamma-1));
if (Pt/Pa) >= PRcrit
    % choked internally; still use fully expanded exit velocity relation
end
Me2 = max(0, (2/(gamma-1))*((Pt/Pa)^((gamma-1)/gamma)-1));
Te = Tt/(1 + (gamma-1)/2*Me2);
Vi = sqrt(max(2*cp*(Tt-Te), 0));
V = Cv*Vi;
end

function [phi, Tmax] = cooling_requirement(des, OPR, atm, M0, eff_bonus, opts)
% Low-order cooling requirement using kSt-scaled transfer and film effectiveness
% (see repo refs: 40_turbine_cooling_annotated.pdf).

porN = (des.porosity - 0.002)/(0.020-0.002);
CdN = (opts.Cd - 0.60)/(0.95-0.60);
angN = sin(deg2rad(opts.alpha));
tN = (des.tTBC - 100)/(400-100);

eta_film = 0.35 + 0.25*porN + 0.15*CdN + 0.10*angN + 0.15*tN + eff_bonus;
eta_film = clamp(eta_film*des.kSt, 0.20, 0.92);

Tg = des.Tt4 - 120*(atm.h_m/11000) + 20*M0;
hScale = des.kSt * (1 + 0.25*((OPR-20)/35) + 0.05*M0);
hScale = max(hScale, 0.4);

phi_req = (Tg - opts.T_allow - 0.35*(des.tTBC - 100)) ...
    /(opts.cool_phi_den*(0.6 + 0.4*eta_film)*hScale);
phi = clamp(phi_req + opts.cool_phi_offset, opts.phi_bounds(1), opts.phi_bounds(2));

Tmax = 980 + 0.55*(des.Tt4 - 1450) ...
    - opts.cool_tmetal_gain*phi*(0.6 + 0.4*eta_film)*hScale ...
    - 0.35*(des.tTBC - 100);
Tmax = clamp(Tmax, opts.Tmetal_floor, opts.Tmetal_ceil);
end

function OPR = OPR_from_des(des)
OPR = des.LPR * des.HPR;
end

function tsfc_out = convert_tsfc_units(tsfc_si, unit_name)
% Convert from SI kg/(N*s) to requested reporting units.
if strcmpi(unit_name, "lbm_per_hr_lbf")
    % 1 kg/(N*s) = 3600*2.20462262185/0.22480894387 lbm/(hr*lbf)
    tsfc_out = tsfc_si * 35303.94013;
else
    tsfc_out = tsfc_si;
end
end

function y = soft_floor(x, x_min, sharpness)
% Smooth floor to avoid hard clipping/flat sensitivity near infeasible states.
y = x_min + (1/sharpness)*log(1 + exp(sharpness*(x - x_min)));
end

function z = clamp(x, lo, hi)
z = min(max(x, lo), hi);
end
