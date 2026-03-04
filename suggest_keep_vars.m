function keep = suggest_keep_vars(mu_star, sigma, varNames, mode)
%SUGGEST_KEEP_VARS Variable selection helper.
% mode options:
%   "hybrid"   : legacy heuristic using top mu* and above-median sigma
%   "pareto80" : keep smallest set reaching 80% cumulative mu* for each
%                response; union the selected sets.

if nargin < 4 || strlength(string(mode)) == 0
    mode = "hybrid";
end
mode = lower(string(mode));

k = numel(varNames);
keep = false(k,1);

switch mode
    case "pareto80"
        nResp = size(mu_star,2);
        for q = 1:nResp
            w = mu_star(:,q);
            w(~isfinite(w)) = 0;
            [ws, idx] = sort(w, 'descend');
            total = sum(ws);
            if total <= 0
                continue;
            end
            cfrac = cumsum(ws) / total;
            nKeep = find(cfrac >= 0.80, 1, 'first');
            if isempty(nKeep); nKeep = k; end
            keep(idx(1:nKeep)) = true;
        end

    otherwise % "hybrid"
        mu_tsfc = mu_star(:,1);
        mu_mt   = mu_star(:,min(3,size(mu_star,2)));

        [~, idx1] = sort(mu_tsfc,'descend');
        [~, idx2] = sort(mu_mt,'descend');

        topN = min(8,k);
        keep(idx1(1:topN)) = true;
        keep(idx2(1:topN)) = true;

        sig_tsfc = sigma(:,1);
        sig_mt   = sigma(:,min(3,size(sigma,2)));
        keep(sig_tsfc > median(sig_tsfc,'omitnan')) = true;
        keep(sig_mt   > median(sig_mt,'omitnan'))   = true;
end
end
