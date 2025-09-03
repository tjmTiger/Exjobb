function test_hubs(tag, cluster_size, start_nodes, colors, graph_name, options)
%TEST_HUBS Summary of this function goes here
%   Detailed explanation goes here
arguments
    tag
    cluster_size
    start_nodes
    colors
    graph_name
    options.n_graphs {mustBeNumeric} = 100
    options.fract_targ_cluster {mustBeNumeric} = 0.1
    options.fract_dist_cluster {mustBeNumeric} = 0.1
end

% get network from data file
load formated_data.mat;
tags = keys(formated_data);
val = values(formated_data);
disp(tags{tag})
G = val{tag}{1}{1};
graph_size = numnodes(G);
% hub_ranks = centrality(G,'hubs');
% G.Nodes.Hubs = hub_ranks;

figure();
fig = plot(G,'Layout','subspace');
fig.NodeColor = '#757575';
fig.EdgeColor = '#a1a1a1';
fig.MarkerSize = 5;
% fig.NodeLabel = 1:numnodes(G);


% [a, b] = find( distances(G) == max(max(distances(G))) );
% disp([a,b])

hubs = [];
for j = 1:length(start_nodes)
    color = colors(j);
    hub = start_nodes(j);
    disp("Hubs left: " + (length(start_nodes)-j))
    while length(hub) < cluster_size
        %temp_hub = hub;
        for i = hub
            p = predecessors(G,i);
            s = successors(G,i);
            highlight(fig,p,'NodeColor',color)
            highlight(fig,s,'NodeColor',color)
            highlight(fig,i,p,'EdgeColor',color)
            highlight(fig,i,s,'EdgeColor',color)
            hub = unique([hub, p',s']);
            if length(hub) >= cluster_size
                break
            end
        end
    end
    hubs = [hubs; hub(1,1:cluster_size)];
end

annotation('textbox',[.78 .9 0 0],'String','\bullet hub1','FitBoxToText','on','Color','r'); % dim: [x y w h]
annotation('textbox',[.78 .835 0 0],'String','\bullet hub2','FitBoxToText','on','Color','b');
annotation('textbox',[.78 .77 0 0],'String','\bullet hub3','FitBoxToText','on','Color','g');
annotation('textbox',[.78 .705 0 0],'String','\bullet hub4','FitBoxToText','on','Color','k');

fontsize(12,"points")
position = get(gcf, 'Position');
position = [100, 100, 600, 600];
saveas(gcf, "figures_new/" + graph_name + " network hubs_graph.fig")

% shortest path
% shortest_paths = [];
% for hub_a = 1:numel(hubs(:,1))
% for hub_b = setdiff(1:numel(hubs(:,1)),hub_a)
%     path_len = 100;
%     for node_a = hubs(hub_a,:)
%     for node_b = hubs(hub_b,:)
%         path_len_temp = length(shortestpath(G,node_a, node_b));
%         if path_len_temp < path_len
%             path_len = path_len_temp;
%         end
%     end
%     end
%     shortest_paths(end+1) = path_len;
% end
% end
% % shortest_paths (Technological 140):   4 5 4 5 5 5
% % shortest_paths (Electrical 600):      2 2 2 2 3 3 2 3 3 2 3 3
% % isempty(intersect(hubs(1,:), hubs(3,:)))

%% solve

figure();
fract_targ = cluster_size*options.fract_targ_cluster/graph_size; % convert from cluster procentage to total procentage
fract_dist = cluster_size*options.fract_dist_cluster/graph_size;

n_graphs = options.n_graphs;
for t = 1:numel(hubs(:,1))
for d = setdiff(1:numel(hubs(:,1)),t)
    disp("hubs: " + t + "-" + d)
    test2do = n_tests(graph_size, fract_targ, fract_dist, cluster_size, cluster_size);

    results = [];
    results_time = [];
    results_trivial = [];
    
    for i = 1:n_graphs
        disp("Left: " + (n_graphs-i) + ", Size: " + graph_size)
        disp("test2do: " + test2do)
        results_i = zeros(1, test2do);
        results_time_i = zeros(1, test2do);
        results_trivial_i = zeros(1, test2do);
        parfor j = 1:test2do
            [results_i(j), results_time_i(j), results_trivial_i(j)] = decouple(G, fract_targ, fract_dist, 'list_targ', hubs(t,:), 'list_dist', hubs(d,:), 'seed', j);
        end
        results = [results results_i];
        results_time = [results_time results_time_i];
        results_trivial = [results_trivial results_trivial_i];
    end
    graph_name = t + "-" + d;% convertCharsToStrings(tags{tag});
    subplot(1,3,1);
    hold on
    add2boxchart(results, graph_name, "Cost", "Cost [-]", "Hubs")
    hold off
    subplot(1,3,2);
    hold on
    add2boxchart(results_time, graph_name, graph_name + " network" + newline + "Runtime", "Time [s]", "Hubs")
    hold off
    subplot(1,3,3);
    hold on
    add2boxchart(results_trivial, graph_name, "Trivial solutions", "Index [-]", "Hubs")
    hold off
    ylim([0 1])
end
end

fontsize(12,"points")
position = get(gcf, 'Position');
position = [100, 100, 600, 600];
saveas(gcf, "figures_new/" + graph_name + " network hubs.fig")
end


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
