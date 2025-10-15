clear;clc;close all;

interpreter = 'latex';
location = "./figures_other/";
local_figures = struct2table(dir(location + "*.fig")).name;

% for figure_name = local_figures'
%     figure_name = cell2mat(figure_name);
%     fig = openfig(location + figure_name);
% 
%     all_text = findall(fig, "-property", "Interpreter");
%     set(all_text, "Interpreter", interpreter)
% 
%     all_text = findall(fig, "-property", "TickLabelInterpreter");
%     set(all_text, "TickLabelInterpreter", interpreter)
% 
%     all_text = findall(fig, "-property", "Fontsize");
%     set(all_text, "Fontsize", 12)
% 
%     savefig(fig, location + figure_name)
%     % print(fig, location +  erase(figure_name,".fig") + ".eps", "-depsc")
%     close all;
% end

% Save all figures in a dir to .eps to specified dir
location2 = "./figures_other/";
for figure_name = local_figures'
    figure_name = cell2mat(figure_name);
    fig = openfig(location + figure_name);

    % savefig(fig, location + figure_name)
    print(fig, location2 +  erase(figure_name,".fig") + ".eps", "-depsc")
    close all;
end

% % for saving legend
% for figure_name = local_figures'
%     figure_name = cell2mat(figure_name);
%     if string(figure_name) == "pareto_ScaleFree_legend.fig"
%         fig = openfig(location + figure_name);
%         set(fig,'Position',[0,0,1024,1024]);
%         legend_handle = fig.Children(2);
%         set(fig,'Position',(get(legend_handle,'Position').*[0, 0, 1, 1].*get(fig,'Position')));
%         set(legend_handle,'Position',[0,0,1,1]);
%         set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);
%     end
% end