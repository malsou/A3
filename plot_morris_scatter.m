function plot_morris_scatter(mu_star, sigma, varNames, ttl, outPng, xLabel)
if nargin < 6 || strlength(string(xLabel)) == 0
    xLabel = '\\mu^*';
end

plotScale = 1;
if max(abs(mu_star)) < 1e-8
    plotScale = 1e6;
    xLabel = string(xLabel) + ' (x10^{-6})';
end
xv = mu_star * plotScale;

fig = figure('Visible','off');
scatter(xv, sigma, 60, 'filled'); grid on;
xlabel(xLabel); ylabel('\\sigma'); title(ttl);

for i = 1:numel(varNames)
    text(xv(i), sigma(i), " " + varNames(i), 'FontSize', 9);
end

set(fig,'Position',[100 100 900 600]);
exportgraphics(fig, outPng, 'Resolution', 200);
close(fig);
end
