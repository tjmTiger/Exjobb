clear;
close all;
clc;

n_graphs = 100; % sample size of graphs
%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

% n_graphs = 100; % number of graphs
for graph_name = ["Erdos Renyi", "Watts Strogatz", "Scale Free"]
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    fract_T_D = 0.05:0.05:0.3;
    for j = fract_T_D
        n = 100;
        fract_targ = j;
        fract_dist = j;
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs
            if graph_name == "Erdos Renyi"
                p = 1.5*(log10(n)/n);
                G = ER_Graph(n, p, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Erdos: " + i)
            elseif graph_name == "Watts Strogatz"
                k = 2;
                beta = 0.2;
                G = WattsStrogatz(n, k, beta, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Strogatz: " + i)
            else
                alpha = 0.2;
                beta = 0.2;
                gamma = 0.6;
                G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("SFG: " + i)
            end
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end

    figure();
    subplot(1,3,1);
    hold on;

    for r = 1:size(results_all, 1)
        add2boxchart(results_all(r,:), string(fract_T_D(r)), "Cost", "Cost [-]", "Fractions")
    end

    hold off;

    subplot(1,3,2);
    hold on;

    for r = 1:size(results_time_all, 1)
        add2boxchart(results_time_all(r,:), string(fract_T_D(r)), graph_name + newline + "Runtime", "Time [s]", "Fractions")
    end

    hold off;

    subplot(1,3,3);
    hold on;

    for r = 1:size(results_trivial_all, 1)
        add2boxchart(results_trivial_all(r,:), string(fract_T_D(r)), "Trivial solutions", "Index [-]", "Fractions")
    end
    ylim([0 1])
    hold off;
end

%-----------------------------------------------%
%                                               %
%                   Graph Size                  %
%                                               %
%-----------------------------------------------%

% n_graphs = 100; % number of graphs
fract_targ = 0.1;
fract_dist = 0.1;

for graph_name = ["Erdos Renyi", "Watts Strogatz", "Scale Free"]
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    graph_size = 30:30:180;
    for j = graph_size
        n = j;
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs
            if graph_name == "Erdos Renyi"
                p = 0.03; % 1.5*(log10(n)/n);
                G = ER_Graph(n, p, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Erdos: " + i)
            elseif graph_name == "Watts Strogatz"
                k = 2;
                beta = 0.2;
                G = WattsStrogatz(n, k, beta, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Strogatz: " + i)
            else
                alpha = 0.2;
                beta = 0.2;
                gamma = 0.6;
                G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("SFG: " + i)
            end
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end

    figure();
    subplot(1,3,1);
    hold on;

    for r = 1:size(results_all, 1)
        add2boxchart(results_all(r,:), string(graph_size(r)), "Cost", "Cost [-]", "Size")
    end

    hold off;

    subplot(1,3,2);
    hold on;

    for r = 1:size(results_time_all, 1)
        add2boxchart(results_time_all(r,:), string(graph_size(r)), graph_name + newline + "Runtime", "Time [s]", "Size")
    end

    hold off;

    subplot(1,3,3);
    hold on;

    for r = 1:size(results_trivial_all, 1)
        add2boxchart(results_trivial_all(r,:), string(graph_size(r)), "Trivial solutions", "Index [-]", "Size")
    end
    ylim([0 1])
    hold off;
end

%-----------------------------------------------%
%                                               %
%                 Connectivity                  %
%                                               %
%-----------------------------------------------%

% n_graphs = 100; % number of graphs
fract_targ = 0.1;
fract_dist = 0.1;
n = 100;
error = 0;
for graph_name = ["Erdos Renyi", "Watts Strogatz", "Scale Free"]
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    node_degree = 2:2:12;
    for j = node_degree
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs
            if graph_name == "Erdos Renyi"
                p = 2.2*j/(n-1); % 1.5*(log10(n)/n);
                G = ER_Graph(n, p, 1, i);
                error = error + (numedges(G)/numnodes(G)-j);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Erdos: " + i)
            elseif graph_name == "Watts Strogatz"
                k = j;
                beta = 0.2;
                G = WattsStrogatz(n, k, beta, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Strogatz: " + i)
            else
                alpha = 0.2;
                beta = 0.2;
                gamma = 0.6;
                G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("SFG: " + i)
            end
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end

    figure();
    subplot(1,3,1);
    hold on;

    for r = 1:size(results_all, 1)
        add2boxchart(results_all(r,:), string(node_degree(r)), "Cost", "Cost [-]", "Average degree")
    end

    hold off;

    subplot(1,3,2);
    hold on;

    for r = 1:size(results_time_all, 1)
        add2boxchart(results_time_all(r,:), string(node_degree(r)), graph_name + newline + "Runtime", "Time [s]", "Average degree")
    end

    hold off;

    subplot(1,3,3);
    hold on;

    for r = 1:size(results_trivial_all, 1)
        add2boxchart(results_trivial_all(r,:), string(node_degree(r)), "Trivial solutions", "Index [-]", "Average degree")
    end
    ylim([0 1])
    hold off;
end

disp("Mean node degree error in Erdos Renyi: " + error/(n_graphs*numel(node_degree))) % -0.0096

%-----------------------------------------------%
%                                               %
%                  Extra tests                  %
%                                               %
%-----------------------------------------------%
%
%-----------------------------------------------%
%                                               %
%                 Watts Strogatz                %
%                                               %
%-----------------------------------------------%
%
% fract_targ = 0.1;
% fract_dist = 0.1;
% 
% results = zeros(1, n_graphs);
% for i = 1:n_graphs
%     G = WattsStrogatz(50, 2, 0.2, 1);
%     results(i) = decouple(G, fract_targ, fract_dist);
% end
% add2boxchart(results, "Watts Strogatz")
%
%-----------------------------------------------%
%                                               %
%                   Scale Free                  %
%                                               %
%-----------------------------------------------%
% fract_targ = 0.1;
% fract_dist = 0.1;
% 
% n_graphs = 100;
% n = 100;
% tot_mean_degree = [];
% results = zeros(1, n_graphs);
% todo = n_graphs*10;
% k = 1;
% figure();
% for beta = [0 0.55 0.70 0.79 0.845 0.871 0.893 0.91 0.923 0.93]
%     for i = 1:n_graphs
%         % beta = 0.1;
%         alpha = (1-beta)/2;
%         gamma = (1-beta)/2;
%         delta_in = 10;
%         delta_out = 10;
%         G = SFG_dir(n, alpha, beta, gamma, delta_in, delta_out, i);
%         [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist); 
%         todo = todo-1;
%         disp("Left: " + todo)
%     end
%     subplot(1,3,1);
%     hold on
%     add2boxchart(results, "$\bar{k}$ = " + k, "Cost", "Cost [-]", "average degree")
%     hold off
%     subplot(1,3,2);
%     hold on
%     add2boxchart(results_time, "$\bar{k}$ = " + k, "Runtime", "Time [s]", "average degree")
%     hold off
%     subplot(1,3,3);
%     hold on
%     add2boxchart(results_trivial, "$\bar{k}$ = " + k, "Trivial solutions", "Index [-]", "average degree")
%     hold off
%     k = k+1;
% end
% add2boxchart(results, "k =" + k_average, "SFG")
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
%     results = zeros(1, n_graphs);
%     results_time = zeros(1, n_graphs);
%     results_trivial = zeros(1, n_graphs);
% 
%     for i = 1:length(val{tag})
%         G = val{tag}{i}{1};
%         disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
%         [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
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
%     graph_name = convertCharsToStrings(tags{tag});
%     hold on
%     add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
%     hold off
% end

% load formated_data.mat;
% 
% tags = keys(formated_data);
% val = values(formated_data);
% % figure();
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
%     end
%     % graph_name = convertCharsToStrings(tags{tag});
%     % hold on
%     % add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
%     % hold off
% end
%
%
% load formated_data.mat;
% 
% tags = keys(formated_data);
% val = values(formated_data);
% for tag = 1:length(tags)
%     disp(tags{tag})
%     n_graphs = length(val{tag});
% 
%     graph_name = convertCharsToStrings(tags{tag});
%     x = ceil(sqrt(n_graphs));
%     figure();
%     for i = 1:length(val{tag})
%         G = val{tag}{i}{1};
%         disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
%         subplot(x,x,i)
%         plot(G)
%     end
%     sgtitle(graph_name)
%     % hold on
%     % add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
%     % hold off
% end
%
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