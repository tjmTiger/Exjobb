clear;clc;close all;

interpreter = 'latex';
location = "./figures_dynamical_feedback_2/";
local_figures = struct2table(dir(location + "*.fig")).name;

for figure_name = local_figures'
    figure_name = cell2mat(figure_name);
    fig = openfig(location + figure_name);

    all_text = findall(fig, "-property", "Interpreter");
    set(all_text, "Interpreter", interpreter)

    all_text = findall(fig, "-property", "TickLabelInterpreter");
    set(all_text, "TickLabelInterpreter", interpreter)

    all_text = findall(fig, "-property", "Fontsize");
    set(all_text, "Fontsize", 12)

    savefig(fig, location + figure_name)
    print(fig, location +  erase(figure_name,".fig") + ".eps", "-depsc")
    close all;
end

% for figure_name = local_figures'
%     figure_name = cell2mat(figure_name);
%     fig = openfig(location + figure_name);
% 
%     % savefig(fig, location + figure_name)
%     print(fig, location +  erase(figure_name,".fig") + ".eps", "-depsc")
%     close all;
% end