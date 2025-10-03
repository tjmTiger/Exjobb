clear;
close all;
clc;
%-----------------------------------------------%
%                                               %
%                  Real Networks                %
%                                               %
%-----------------------------------------------%
% 
% fract_targ = 0.1;
% fract_dist = 0.1;
% 
% load formated_data.mat;
% 
% tags = keys(formated_data);
% val = values(formated_data);
% figure();
% results_all = [];
% results_time_all = [];
% results_trivial_all = [];
% for tag = 1:length(tags)
%     disp(tags{tag})
% 
%     n_graphs = length(val{tag});
%     results = [];
%     results_time = [];
%     results_trivial = [];
% 
%     for i = 1:length(val{tag})
%         G = val{tag}{i}{1};
%         disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
%         test2do = n_tests(size(G.Nodes, 1), fract_targ, fract_dist);
%         disp("test2do: " + test2do)
%         for j = 1:test2do
%             [results(end+1), results_time(end+1), results_trivial(end+1)] = decouple(G, fract_targ, fract_dist);
%         end
%     end
%     graph_name = convertCharsToStrings(tags{tag});
%     subplot(1,3,1);
%     hold on
%     add2boxchart(results, graph_name, "Cost", "Cost [-]", "Graph category")
%     hold off
%     subplot(1,3,2);
%     hold on
%     add2boxchart(results_time, graph_name, "Runtime", "Time [s]", "Graph category")
%     hold off
%     subplot(1,3,3);
%     hold on
%     add2boxchart(results_trivial, graph_name, "Trivial solutions", "Index [-]", "Graph category")
%     hold off
% end
% 
% fontsize(12,"points")
% position = get(gcf, 'Position');
% position = [100, 100, 600, 600];
% saveas(gcf, "figures_new/Real Networks.fig")

% %% Sizes of networks
% load formated_data.mat;
% 
% tags = keys(formated_data);
% val = values(formated_data);
% figure();
% for tag = 1:length(tags)
%     disp(tags{tag})
% 
%     n_graphs = length(val{tag});
%     results = zeros(1, n_graphs);
%     results_time = zeros(1, n_graphs);
%     results_trivial = zeros(1, n_graphs);
% 
%     for i = 1:length(val{tag})
%         G = val{tag}{i}{1};
%         disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
%         results(i) = size(G.Nodes, 1);
%     end
%     graph_name = erase(convertCharsToStrings(tags{tag}), " network");
%     graph_name = replace(graph_name, "communication", "com.");
%     graph_name = replace(graph_name, "communitie", "com.");
%     hold on
%     add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
%     hold off
%     set(gca, 'YScale', 'log')
% end
%
%% plot all graph topologies
load formated_data.mat;

tags = keys(formated_data);
val = values(formated_data);
for tag = 1:length(tags)
    disp(tags{tag})
    n_graphs = length(val{tag});

    graph_name = convertCharsToStrings(tags{tag});
    x = ceil(sqrt(n_graphs));
    % figure();
    for i = 1:length(val{tag})
        G = val{tag}{i}{1};
        disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
        % subplot(x*x,1,i)
        figure();
        plot(G)
        print(gcf, "figures_new/" + erase(graph_name," ") + "_" + num2str(i) + ".eps", "-depsc")
        close all;
    end
    sgtitle(graph_name)
end

% %%
% %-----------------------------------------------%
% %                                               %
% %              Technological network            %
% %                                               %
% %-----------------------------------------------%
% clear;
% clc;
% multithreading()
% 
% % number of nodes: 677
% test_hubs(11, 140, [340, 80, 359], ['r', 'b', 'g'], "Technological")
% 
% %
% %-----------------------------------------------%
% %                                               %
% %                Electrical network             %
% %                                               %
% %-----------------------------------------------%
% clear;
% clc;
% multithreading()
% 
% % number of nodes: 4941
% test_hubs(12, 500, [4463, 3743, 2539, 97], ['r', 'b', 'g','k'], "Electrical")

%-----------------------------------------------%
%                                               %
%                   Functions                   %
%                                               %
%-----------------------------------------------%

function add2boxchart(results, test_name, title_name, ylabel_name, xlabel_name)
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    boxchart(categorical(1:numel(results), 1:numel(results), repmat(test_name, 1, numel(results))), results)
    if nargin >= 3
        title(title_name)
    end
    if nargin == 4
        ylabel(ylabel_name)
        xlabel('Tests')
    elseif nargin == 5
        ylabel(ylabel_name)
        xlabel(xlabel_name)
    else
        ylabel('Cost [-]')
        xlabel('Tests')
    end
end

function n_tests = n_tests(n, fracT, fracD, hubT, hubD)
    alpha  = 1;
    switch nargin
        case 3
            n_tests = alpha*log(n)*1/(fracT+fracD);
        case 5
            n_tests = alpha*log(hubT+hubD)*(hubT/fracT + hubD/fracD)*(1/n);
        otherwise
            disp('input argument invalid')
    end
    n_tests = min(100, max(10,floor(n_tests)));
end

function [a,b,c] = test()
    a = 1; b = 1; c = 1;
end