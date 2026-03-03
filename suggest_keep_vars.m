function keep = suggest_keep_vars(mu_star, sigma, varNames)
% Simple rule:
% Keep if in top 8 by mu* for TSFC or MT, OR sigma above median for either.

k = numel(varNames);
keep = false(k,1);

mu_tsfc = mu_star(:,1);
mu_mt   = mu_star(:,3);

[~, idx1] = sort(mu_tsfc,'descend');
[~, idx2] = sort(mu_mt,'descend');

topN = min(8,k);
keep(idx1(1:topN)) = true;
keep(idx2(1:topN)) = true;

sig_tsfc = sigma(:,1);
sig_mt   = sigma(:,3);
keep(sig_tsfc > median(sig_tsfc,'omitnan')) = true;
keep(sig_mt   > median(sig_mt,'omitnan'))   = true;
end
