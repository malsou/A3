%% A3: Screening setup + Morris screening (no toolboxes)
% Run this script. Outputs go to ./a3_out/
clear; clc; close all;
rng(7);

outDir = "a3_out";
if ~exist(outDir, "dir"); mkdir(outDir); end

%% ----- Design space (A2/A3) -----
% Categorical variables (handled as cases)
nozzleCases  = ["mixed","separate"];    % x1
coolingCases = ["film","effusion"];     % x2

% Continuous variables (k = 11)
% NOTE: eta_pc and eta_pt are no longer independent DOE inputs.
% They are computed inside evaluate_engine_model_proxy as coupled functions of
% cycle + cooling state, while SMc and kSt are independent knobs.
varNames = ["BPR","FPR","LPR","HPR","Tt4","SMc","kSt","porosity","Cd","alpha","tTBC"];

% Units (for tables)
varUnits = ["-","-","-","-","K","%","-","-","-","deg","um"];

% Bounds (keep yours; I recommend FPR min >= 1.30 for feasibility)
lb = [ 4.0, 1.30, 1.20,  8.0, 1250, 0.05, 0.70, 0.002, 0.60, 20, 100 ];
ub = [12.0, 1.90, 3.00, 30.0, 1900, 0.25, 1.50, 0.020, 0.95, 90, 400 ];

% Derived feasibility constraint for OPR = LPR*HPR
OPR_min = 20;
OPR_max = 55;

%% ----- Morris design settings -----
k = numel(varNames);
p = 6;         % grid levels
r = 8;         % trajectories -> N = r*(k+1) runs per case (96)
delta = p/(2*(p-1));  % standard Morris step

fprintf("Morris settings: k=%d, p=%d, r=%d, delta=%.3f, N=%d per case\n", ...
    k,p,r,delta,r*(k+1));

%% ----- Choose screening strategy for categorical variables -----
% Option A: baseline case Morris + small categorical comparison
DO_BASELINE_ONLY = true;
baselineNozzle = "mixed";
baselineCooling = "film";

% Small categorical comparison sample size (optional)
N_cat_compare = 10;

%% ----- Generate Morris design in normalized space -----
% Keep only trajectories that satisfy OPR constraints at all trajectory points
% to preserve valid one-factor-at-a-time Morris steps.
denorm = @(Xn) denormalize_minmax(Xn, lb, ub);
iL = find(varNames=="LPR",1);
iH = find(varNames=="HPR",1);

Xn_all = zeros(r*(k+1), k);
trajAccepted = 0;
trajAttempts = 0;
maxTrajAttempts = 5000;

while trajAccepted < r && trajAttempts < maxTrajAttempts
    trajAttempts = trajAttempts + 1;
    Xt = morris_design(k,p,1);               % one trajectory: (k+1) x k
    Xraw_t = denorm(Xt);
    opr_t = Xraw_t(:,iL) .* Xraw_t(:,iH);

    if all(opr_t >= OPR_min & opr_t <= OPR_max)
        rowStart = trajAccepted*(k+1) + 1;
        Xn_all(rowStart:rowStart+k, :) = Xt;
        trajAccepted = trajAccepted + 1;
    end
end

if trajAccepted < r
    error("Could not generate %d feasible Morris trajectories in %d attempts.", r, maxTrajAttempts);
end
fprintf("Accepted %d/%d feasible Morris trajectories after %d attempts.\n", r, r, trajAttempts);

%% ----- Run screening -----
results = struct;

caseList = {};
if DO_BASELINE_ONLY
    caseList = {baselineNozzle, baselineCooling};
else
    for a = 1:numel(nozzleCases)
        for b = 1:numel(coolingCases)
            caseList(end+1,:) = {nozzleCases(a), coolingCases(b)}; 
        end
    end
end

for ci = 1:size(caseList,1)
    nozzle = string(caseList{ci,1});
    cooling = string(caseList{ci,2});
    caseTag = nozzle + "_" + cooling;

    fprintf("\n=== Running Morris screening case: %s ===\n", caseTag);

    Xraw = denorm(Xn_all);

    % Evaluate
    N = size(Xraw,1);
    Y = zeros(N,4); % [TSFC_cruise, phi_crit, MT_worst, MF_worst]
    runtime = zeros(N,1);

    for i = 1:N
        x = row_to_struct(Xraw(i,:), varNames);
        x.nozzle  = nozzle;
        x.cooling = cooling;

        t0 = tic;
        y = evaluate_engine_model_proxy(x);
        runtime(i) = toc(t0);

        Y(i,:) = [y.TSFC_cruise, y.phi_crit, y.MT_worst, y.MF_worst];
    end

    % Save dataset CSV
    dataVarNames = [varNames, "TSFC_cruise","phi_crit","MT_worst","MF_worst","runtime_s"];
    T = array2table([Xraw, Y, runtime], "VariableNames", cellstr(dataVarNames));
    writetable(T, fullfile(outDir, "screening_dataset_" + caseTag + ".csv"));

    % Morris analysis (elementary effects)
    mor = morris_analyze(Xn_all, Y, k, r, delta);

    % Plots
    plot_morris_bars(mor.mu_star(:,1), varNames, "Morris \\mu^* (TSFC)", fullfile(outDir, "morris_mu_TSFC_" + caseTag + ".png"));
    plot_morris_bars(mor.mu_star(:,3), varNames, "Morris \\mu^* (MT margin)", fullfile(outDir, "morris_mu_MT_" + caseTag + ".png"));

    plot_morris_scatter(mor.mu_star(:,1), mor.sigma(:,1), varNames, ...
        "Morris: TSFC ( \\mu^* vs \\sigma )", fullfile(outDir, "morris_scatter_TSFC_" + caseTag + ".png"));
    plot_morris_scatter(mor.mu_star(:,3), mor.sigma(:,3), varNames, ...
        "Morris: MT margin ( \\mu^* vs \\sigma )", fullfile(outDir, "morris_scatter_MT_" + caseTag + ".png"));

    % Variable keep/drop suggestion (simple rule)
    keep = suggest_keep_vars(mor.mu_star, mor.sigma, varNames);

    % Save ranking table
    R = table(varNames', mor.mu_star(:,1), mor.sigma(:,1), mor.mu_star(:,3), mor.sigma(:,3), keep, ...
        'VariableNames', ["Variable","mu_TSFC","sigma_TSFC","mu_MT","sigma_MT","Keep"]);
    writetable(R, fullfile(outDir, "screening_rank_" + caseTag + ".csv"));

    % Summary print
    fprintf("Avg runtime per eval: %.4f s (N=%d)\n", mean(runtime), N);
    disp(R);

    results.(caseTag).dataset = T;
    results.(caseTag).morris = mor;
    results.(caseTag).rank = R;
    results.(caseTag).avg_runtime = mean(runtime);
end

%% ----- Optional: categorical comparison (small LHS sample replicated across categories) -----
if DO_BASELINE_ONLY
    fprintf("\n=== Optional categorical comparison (N=%d) ===\n", N_cat_compare);
    Xn_cat = lhs_simple(N_cat_compare, k); % simple LHS in [0,1]
    Xraw_cat = denorm(Xn_cat);
    Xraw_cat = enforce_opr_bounds(Xraw_cat, varNames, lb, ub, OPR_min, OPR_max);

    nCatRows = numel(nozzleCases) * numel(coolingCases) * N_cat_compare;
    nozzleCol = strings(nCatRows,1);
    coolingCol = strings(nCatRows,1);
    Xcat = zeros(nCatRows, k);
    Ycat = zeros(nCatRows, 4);

    rr = 0;
    for a = 1:numel(nozzleCases)
        for b = 1:numel(coolingCases)
            nozzle = nozzleCases(a);
            cooling = coolingCases(b);
            for i = 1:N_cat_compare
                rr = rr + 1;
                x = row_to_struct(Xraw_cat(i,:), varNames);
                x.nozzle = nozzle;
                x.cooling = cooling;
                y = evaluate_engine_model_proxy(x);

                nozzleCol(rr) = nozzle;
                coolingCol(rr) = cooling;
                Xcat(rr,:) = Xraw_cat(i,:);
                Ycat(rr,:) = [y.TSFC_cruise, y.phi_crit, y.MT_worst, y.MF_worst];
            end
        end
    end

    Tcat = table(nozzleCol, coolingCol, 'VariableNames', {'nozzle','cooling'});
    Tcont = array2table(Xcat, 'VariableNames', cellstr(varNames));
    Tout = array2table(Ycat, 'VariableNames', {'TSFC_cruise','phi_crit','MT_worst','MF_worst'});
    Tcat = [Tcat, Tcont, Tout];
    writetable(Tcat, fullfile(outDir, "categorical_compare.csv"));
end

%% ----- Export a small normalization example table (for report) -----
exampleRow = 1;
Xraw_ex = denorm(Xn_all(exampleRow,:));
Xn_ex = normalize_minmax(Xraw_ex, lb, ub);

export_norm_example_tex(varNames, varUnits, Xraw_ex, lb, ub, Xn_ex, fullfile(outDir, "norm_example_table.tex"));

fprintf("\nDone. Outputs in: %s\n", outDir);
