%% A3: Screening setup + Morris screening (no toolboxes)
% Categorical factors are treated as distinct scenarios.
% Morris is run only on continuous variables within each scenario.
clear; clc; close all;
rng(7);

outDir = "a3_out";
if ~exist(outDir, "dir"); mkdir(outDir); end

%% ----- Design space (continuous + categorical scenarios) -----
nozzleCases  = ["mixed","separate"];
coolingCases = ["film","effusion"];

% Continuous variables only (Morris dimensions)
varNames = ["BPR","FPR","LPR","HPR","Tt4","SMc","kSt","porosity","Cd","alpha","tTBC"];
varUnits = ["-","-","-","-","K","%","-","-","-","deg","um"];

lb = [ 4.0, 1.30, 1.20,  8.0, 1250, 0.05, 0.70, 0.002, 0.60, 20, 100 ];
ub = [12.0, 1.90, 3.00, 30.0, 1900, 0.25, 1.50, 0.020, 0.95, 90, 400 ];

OPR_min = 20;
OPR_max = 55;

%% ----- Morris design settings -----
k = numel(varNames);
p = 6;
r = 8;
delta = p/(2*(p-1));

fprintf("Morris settings (continuous only): k=%d, p=%d, r=%d, delta=%.3f, N=%d per scenario\n", ...
    k,p,r,delta,r*(k+1));

%% ----- Build one feasible continuous Morris design -----
denorm = @(Xn) denormalize_minmax(Xn, lb, ub);
iL = find(varNames=="LPR",1);
iH = find(varNames=="HPR",1);

Xn_all = zeros(r*(k+1), k);
trajAccepted = 0;
trajAttempts = 0;
maxTrajAttempts = 5000;

while trajAccepted < r && trajAttempts < maxTrajAttempts
    trajAttempts = trajAttempts + 1;
    Xt = morris_design(k,p,1);
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
fprintf("Accepted %d/%d feasible trajectories after %d attempts.\n", r, r, trajAttempts);

Xraw = denorm(Xn_all);

%% ----- Scenario list for categorical factors -----
scenarios = struct('nozzle',{},'cooling',{},'tag',{});
for a = 1:numel(nozzleCases)
    for b = 1:numel(coolingCases)
        tag = nozzleCases(a) + "_" + coolingCases(b);
        scenarios(end+1) = struct('nozzle',nozzleCases(a),'cooling',coolingCases(b),'tag',tag); %#ok<SAGROW>
    end
end

%% ----- Run Morris per categorical scenario -----
results = struct;
summaryRows = table();

for si = 1:numel(scenarios)
    nozzle = scenarios(si).nozzle;
    cooling = scenarios(si).cooling;
    caseTag = scenarios(si).tag;

    scenDir = fullfile(outDir, caseTag);
    if ~exist(scenDir, "dir"); mkdir(scenDir); end

    fprintf("\n=== Scenario: %s (nozzle=%s, cooling=%s) ===\n", caseTag, nozzle, cooling);

    N = size(Xraw,1);
    Y = zeros(N,4); % [TSFC_cruise, phi_crit, MT_worst, MF_worst]
    runtime = zeros(N,1);

    for i = 1:N
        x = row_to_struct(Xraw(i,:), varNames);
        x.nozzle = nozzle;
        x.cooling = cooling;

        t0 = tic;
        y = evaluate_engine_model_proxy(x);
        runtime(i) = toc(t0);

        % Keep responses continuous (no pass/fail binning).
        Y(i,:) = [y.TSFC_cruise, y.phi_crit, y.MT_worst, y.MF_worst];
    end

    dataVarNames = [varNames, "TSFC_cruise","phi_crit","MT_worst","MF_worst","runtime_s"];
    T = array2table([Xraw, Y, runtime], "VariableNames", cellstr(dataVarNames));
    writetable(T, fullfile(scenDir, "screening_dataset_" + caseTag + ".csv"));

    % Morris: mu (signed), mu_star, sigma2 (= var(di,1)), sigma.
    mor = morris_analyze(Xn_all, Y, k, r, delta);

    plot_morris_bars(mor.mu_star(:,1), varNames, "Morris \\mu^* (TSFC)", fullfile(scenDir, "morris_mu_TSFC_" + caseTag + ".png"));
    plot_morris_bars(mor.mu_star(:,3), varNames, "Morris \\mu^* (MT margin)", fullfile(scenDir, "morris_mu_MT_" + caseTag + ".png"));

    plot_morris_scatter(mor.mu_star(:,1), mor.sigma(:,1), varNames, ...
        "Morris: TSFC ( \\mu^* vs \\sigma )", fullfile(scenDir, "morris_scatter_TSFC_" + caseTag + ".png"));
    plot_morris_scatter(mor.mu_star(:,3), mor.sigma(:,3), varNames, ...
        "Morris: MT margin ( \\mu^* vs \\sigma )", fullfile(scenDir, "morris_scatter_MT_" + caseTag + ".png"));

    keep = suggest_keep_vars(mor.mu_star, mor.sigma, varNames, "pareto80");

    R = table(varNames', mor.mu(:,1), mor.mu_star(:,1), mor.sigma2(:,1), mor.sigma(:,1), ...
              mor.mu(:,3), mor.mu_star(:,3), mor.sigma2(:,3), mor.sigma(:,3), keep, ...
        'VariableNames', ["Variable","mu_TSFC","mu_star_TSFC","sigma2_TSFC","sigma_TSFC", ...
                          "mu_MT","mu_star_MT","sigma2_MT","sigma_MT","Keep"]);
    writetable(R, fullfile(scenDir, "screening_rank_" + caseTag + ".csv"));

    % Aggregate summary row: top variables and runtime
    [~, idxTSFC] = max(mor.mu_star(:,1));
    [~, idxMT] = max(mor.mu_star(:,3));
    summaryRows = [summaryRows; table(caseTag, string(varNames(idxTSFC)), string(varNames(idxMT)), ...
                    mean(runtime), sum(keep), ...
                    'VariableNames', {"scenario","top_TSFC_mu_star","top_MT_mu_star","avg_runtime_s","n_keep_pareto80"})]; %#ok<AGROW>

    fprintf("Avg runtime per eval: %.4f s (N=%d)\n", mean(runtime), N);
    disp(R);

    results.(caseTag).dataset = T;
    results.(caseTag).morris = mor;
    results.(caseTag).rank = R;
    results.(caseTag).avg_runtime = mean(runtime);
end

% Aggregate comparison across categorical scenarios
writetable(summaryRows, fullfile(outDir, "screening_rank_summary_across_scenarios.csv"));

%% ----- Export normalization example table (~5 rows from same samples) -----
nEx = min(5, size(Xn_all,1));
Xraw_ex = Xraw(1:nEx,:);
Xn_ex = normalize_minmax(Xraw_ex, lb, ub);
export_norm_example_tex(varNames, varUnits, Xraw_ex, lb, ub, Xn_ex, fullfile(outDir, "norm_example_table.tex"));

fprintf("\nDone. Outputs in: %s\n", outDir);

function X = denormalize_minmax(Xn, lb, ub)
X = lb + Xn .* (ub - lb);
end

function Xn = normalize_minmax(X, lb, ub)
Xn = (X - lb) ./ (ub - lb);
end

function s = row_to_struct(row, varNames)
s = struct;
for j = 1:numel(varNames)
    s.(varNames(j)) = row(j);
end
end

function X = lhs_simple(n, k)
% Latin hypercube in [0,1] (basic)
X = zeros(n,k);
for j = 1:k
    edges = linspace(0,1,n+1);
    u = rand(n,1);
    X(:,j) = edges(1:n)' + u.*(edges(2:n+1)'-edges(1:n)');
    X(:,j) = X(randperm(n),j);
end
end

function X = enforce_opr_bounds(X, varNames, lb, ub, OPR_min, OPR_max)
% Resample LPR/HPR for rows violating OPR bounds.
% Throws an error if a row cannot be repaired within maxTries.

iL = find(varNames=="LPR",1);
iH = find(varNames=="HPR",1);
if isempty(iL) || isempty(iH)
    error("LPR/HPR not found");
end

maxTries = 200;
for i = 1:size(X,1)
    opr = X(i,iL) * X(i,iH);
    tries = 0;
    while (opr < OPR_min || opr > OPR_max) && tries < maxTries
        % resample LPR/HPR uniformly within bounds
        X(i,iL) = lb(iL) + rand*(ub(iL)-lb(iL));
        X(i,iH) = lb(iH) + rand*(ub(iH)-lb(iH));
        opr = X(i,iL) * X(i,iH);
        tries = tries + 1;
    end

    if opr < OPR_min || opr > OPR_max
        error("Failed to enforce OPR bounds for row %d after %d retries (OPR=%.4f).", i, maxTries, opr);
    end
end
end

function export_norm_example_tex(varNames, varUnits, xRaw, lb, ub, xN, outTex)
fid = fopen(outTex,'w');
fprintf(fid,"%% Auto-generated normalization example table\n");
fprintf(fid,"\\begin{table}[H]\n\\centering\n");
fprintf(fid,"\\caption{Example pre- and post-normalization values (min--max).}\\label{tab:norm_example}\n");
fprintf(fid,"\\renewcommand{\\arraystretch}{1.15}\n");
fprintf(fid,"\\begin{tabular}{lccccc}\n\\toprule\n");
fprintf(fid,"\\textbf{Var} & \\textbf{Units} & \\textbf{Raw} & \\textbf{Min} & \\textbf{Max} & \\textbf{Norm}\\\\\n\\midrule\n");

for i = 1:numel(varNames)
    fprintf(fid,"%s & %s & %.4g & %.4g & %.4g & %.4g\\\\\n", ...
        varNames(i), varUnits(i), xRaw(i), lb(i), ub(i), xN(i));
end

fprintf(fid,"\\bottomrule\n\\end{tabular}\n\\end{table}\n");
fclose(fid);
end
