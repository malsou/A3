function plot_morris_bars(mu_star, varNames, ttl, outPng)
[vals, idx] = sort(mu_star, 'descend');
names = varNames(idx);

fig = figure('Visible','off');
bar(vals);
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',45);
ylabel('\\mu^*'); title(ttl);
grid on;
set(fig,'Position',[100 100 1100 450]);
exportgraphics(fig, outPng, 'Resolution', 200);
close(fig);
end
