function opts = engine_defaults(overrides)
%ENGINE_DEFAULTS Fixed assumptions for eval_engine.

opts = struct;

% Coupled compressor efficiency model constants
opts.eta_pc_peak = 0.91;
opts.OPR0 = 35;
opts.SMc0 = 0.15;
opts.eta_pc_a = 0.090;
opts.eta_pc_b = 1.80;
opts.eta_pc_bounds = [0.80, 0.92];
opts.SMc = 0.15;

% Turbine efficiency coupling
opts.eta_pt_peak = 0.94;
opts.eta_pt_penalty_c = 1.4;
opts.eta_pt_bounds = [0.82, 0.94];

% Nozzle/cooling categorical modifiers
opts.nozFactorMixed = -0.015;   % mixed nozzle TSFC benefit
opts.coolEffEffusion = 0.08;    % effusion effectiveness bonus

% Cycle proxy constants
opts.ST_base = 240;
opts.ST_Tt4_gain = 120/300;
opts.ST_BPR_gain = -28/4;
opts.ST_OPR_gain = -18/15;
opts.ST_min = 80;
opts.mdot_core_min = 50;

opts.eta_th_bounds = [0.20, 0.55];
opts.eta_p_bounds = [0.40, 0.80];

% Cooling / thermal model constants (in spirit of course cooling notes)
opts.Cd = 0.80;
opts.alpha = 45;                % deg, injection angle proxy
opts.T_allow = 1100;            % K
opts.Tmetal_floor = 750;
opts.Tmetal_ceil = 1400;
opts.phi_bounds = [0.0, 0.08];

% Margin / objective constants
opts.J_weights = [0.75, 0.25];  % [cruise, climb]
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
