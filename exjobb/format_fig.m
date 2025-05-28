clear;
close all;
clc;

graph_names = ["erdos", "sfg", "strogatz"];
test_names = ["fract", "node_degree", "size"];
size = 12;
position = [100, 100, 600, 600]; % x, y, width, height

for graph_name = graph_names
    for test_name = test_names
        figure_name = graph_name + "_" + test_name;

        fig = openfig("figures/" + figure_name + ".fig");
        ax_handles = findall(fig, 'Type', 'Axes');
        ax_handles(1).YLim = [0, 1];
        for i = 1:3
            ax_handles(i).FontSize = size;
            ax_handles(i).Position(1) = ax_handles(i).Position(1) + abs(ax_handles(i).Position(1)*0.1);
            % ax_handles(i).YLabel.Position(1) = ax_handles(i).YLabel.Position(1) + (ax_handles(i).YLabel.Position(1) - ax_handles(i).XLim(1))*0.5;
            % ax_handles(i).YLabel.Position(2) = ax_handles(i).YLabel.Position(2) + (ax_handles(i).YLabel.Position(2) - ax_handles(i).YLim(1))*1.05;
            % ax_handles(i).YLabel.Rotation = 0;
        end

        set(fig, 'Position',  position)
        
        % savefig("figures_new/" + figure_name + ".fig")
        saveas(gcf, "figures_new/" + figure_name + ".jpg")
        close all;
    end
end

fig = openfig("figures/real_networks.fig");
ax_handles = findall(fig, 'Type', 'Axes');
% ax_handles(1).YLim = [-1, 1];
ax_handles(1).FontSize = size;
ax_handles(2).FontSize = size;
% ax_handles(3).FontSize = size;
set(fig, 'Position',  position)

saveas(gcf, "figures_new/real_networks.jpg")
close all;

fig = openfig("figures/real_networks_sizes.fig");
ax_handles = findall(fig, 'Type', 'Axes');
% ax_handles(1).YLim = [-1, 1];
ax_handles(1).FontSize = size-1;
% ax_handles(2).FontSize = size;
% ax_handles(3).FontSize = size;
set(fig, 'Position',  position)

saveas(gcf, "figures_new/real_networks_sizes.jpg")
close all;

fig = openfig("figures/real_networks_trivial_sol.fig");
ax_handles = findall(fig, 'Type', 'Axes');
% ax_handles(1).YLim = [-1, 1];
ax_handles(1).FontSize = size;
ax_handles(2).FontSize = size;
ax_handles(3).FontSize = size;
set(fig, 'Position',  position)

saveas(gcf, "figures_new/real_networks_trivial_sol.jpg")
close all;