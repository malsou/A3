function plot_morris_bars(mu_star, varNames, ttl, outPng)
[vals, idx] = sort(mu_star, 'descend');
names = varNames(idx);

plotScale = 1;
unitTag = "";
if max(abs(vals)) < 1e-8
    plotScale = 1e6;
    unitTag = " (x10^{-6})";
end

fig = figure('Visible','off');
bar(vals*plotScale);
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',45);
ylabel("\\mu^*" + unitTag); title(ttl);
grid on;
set(fig,'Position',[100 100 1100 450]);
exportgraphics(fig, outPng, 'Resolution', 200);
close(fig);
end
