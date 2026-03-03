function plot_morris_scatter(mu_star, sigma, varNames, ttl, outPng)
fig = figure('Visible','off');
scatter(mu_star, sigma, 60, 'filled'); grid on;
xlabel('\\mu^*'); ylabel('\\sigma'); title(ttl);

for i = 1:numel(varNames)
    text(mu_star(i), sigma(i), " " + varNames(i), 'FontSize', 9);
end

set(fig,'Position',[100 100 900 600]);
exportgraphics(fig, outPng, 'Resolution', 200);
close(fig);
end
