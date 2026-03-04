function opts = engine_defaults(overrides)
%ENGINE_DEFAULTS Fixed assumptions and constants for eval_engine.

opts = struct;

% Gas / fuel
opts.gamma_cold = 1.40;
opts.gamma_hot = 1.3333;
opts.cp_cold = 1004.5;          % J/(kg-K)
opts.cp_hot = 1150.0;           % J/(kg-K)
opts.R = 287.05287;             % J/(kg-K)
opts.FHV = 43.1e6;              % J/kg

% Component models
opts.eta_r = 0.99;              % inlet pressure recovery
opts.eta_fan = 0.88;
opts.eta_comp_peak = 0.90;      % baseline for coupled eta_pc
opts.eta_b = 0.995;
opts.dPt_comb = 0.05;           % combustor pressure-loss fraction

% Coupled compressor efficiency model
opts.OPR0 = 35;
opts.SMc0 = 0.15;
opts.eta_pc_a = 0.090;
opts.eta_pc_b = 1.80;
opts.eta_pc_bounds = [0.82, 0.91];
opts.SMc = 0.15;

% Turbine model
opts.eta_t_uncooled = 0.93;
opts.eta_pt_penalty_c = 1.4;    % eta_pt = eta_t_uncooled - c*phi
opts.eta_pt_bounds = [0.84, 0.94];

% Bleeds / cooling baseline fractions (of core flow)
opts.fr_bleed = 0.03;
opts.fr_leak = 0.01;
opts.fr_cool_base = 0.02;

% Nozzle velocity coefficients
opts.Cv19 = 0.985;
opts.Cv9 = 0.983;

% Cooling constants
opts.Cd = 0.80;
opts.alpha = 45;                % deg
opts.T_allow = 1100;            % K
opts.phi_bounds = [0.0, 0.20];
opts.phi_scale = 6000;          % cooling scaling for phi_req and Tmax coupling
opts.Tmetal_floor = 700;
opts.Tmetal_ceil = 1500;
opts.coolEffEffusion = 0.08;
opts.nozFactorMixed = -0.010;

% Objective / penalty
opts.J_weights = [0.75, 0.25];
opts.Tscale = 100;
opts.beta = 1.0;
opts.gamma = 10.0;
opts.P_solver = 0;

if nargin >= 1 && ~isempty(overrides)
    f = fieldnames(overrides);
    for i = 1:numel(f)
        opts.(f{i}) = overrides.(f{i});
    end
end
end
