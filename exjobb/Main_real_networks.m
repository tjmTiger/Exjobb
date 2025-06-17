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
%     graph_name = convertCharsToStrings(tags{tag});
%     hold on
%     add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
%     hold off
%     set(gca, 'YScale', 'log')
% end
%
% %% plot all graph topologies
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
% end

%%
%-----------------------------------------------%
%                                               %
%              Technological network            %
%                                               %
%-----------------------------------------------%
clear;
clc;

load formated_data.mat;
tags = keys(formated_data);
val = values(formated_data);
tag = 11; % Technological network


disp(tags{tag})
G = val{tag}{1}{1};
hub_ranks = centrality(G,'hubs');
% G.Nodes.Hubs = hub_ranks;

figure();
fig = plot(G,'Layout','subspace');
fig.NodeColor = '#757575';
fig.EdgeColor = '#a1a1a1';
fig.MarkerSize = 5;
% fig.NodeLabel = 1:numnodes(G);


% [a, b] = find( distances(G) == max(max(distances(G))) );
% disp([a,b])

cluster_size = 140;
hubs = [];
start_nodes = [340, 80, 359];
colors = ['r', 'b', 'g'];

for j = 1:length(start_nodes)
    color = colors(j);
    hub = start_nodes(j);
    while length(hub) < cluster_size
        temp_hub = hub;
        for i = temp_hub
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

fontsize(12,"points")
position = get(gcf, 'Position');
position = [100, 100, 600, 600];
saveas(gcf, "figures_new/Technological network hubs_graph.fig")

figure();
fract_targ = 0.1;
fract_dist = 0.1;

n_graphs = 100;

for t = 1:numel(hubs(:,1))
for d = setdiff(1:numel(hubs(:,1)),t)
    results = [];
    results_time = [];
    results_trivial = [];
    
    for i = 1:n_graphs
        disp("Left: " + (n_graphs-i) + ", Size: " + size(G.Nodes, 1))
        test2do = n_tests(size(G.Nodes, 1), fract_targ, fract_dist, cluster_size, cluster_size);
        disp("test2do: " + test2do)
        for j = 1:test2do
            [results(end+1), results_time(end+1), results_trivial(end+1)] = decouple(G, fract_targ, fract_dist, hubs(t,:), hubs(d,:));
        end
    end
    graph_name = t + "-" + d;% convertCharsToStrings(tags{tag});
    subplot(1,3,1);
    hold on
    add2boxchart(results, graph_name, "Cost", "Cost [-]", "Hubs")
    hold off
    subplot(1,3,2);
    hold on
    add2boxchart(results_time, graph_name, "Technological network" + newline + "Runtime", "Time [s]", "Hubs")
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
saveas(gcf, "figures_new/Technological network hubs.fig")

%%
%-----------------------------------------------%
%                                               %
%                Electrical network             %
%                                               %
%-----------------------------------------------%
clear;
clc;

load formated_data.mat;
tags = keys(formated_data);
val = values(formated_data);
tag = 12; % Electrical network


disp(tags{tag})
G = val{tag}{1}{1};
hub_ranks = centrality(G,'hubs');
% G.Nodes.Hubs = hub_ranks;

figure();
fig = plot(G,'Layout','subspace');
fig.NodeColor = '#757575';
fig.EdgeColor = '#a1a1a1';
fig.MarkerSize = 5;
% fig.NodeLabel = 1:numnodes(G);


% [a, b] = find( distances(G) == max(max(distances(G))) );
% disp([a,b])

cluster_size = 600;
hubs = [];
start_nodes = [4463, 3743, 2539, 97]; % 4463
colors = ['r', 'b', 'g','k'];

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
saveas(gcf, "figures_new/Electrical network hubs_graph.fig")

%%

figure();
fract_targ = 0.1;
fract_dist = 0.1;

n_graphs = 10;

for t = 1:numel(hubs(:,1))
for d = setdiff(1:numel(hubs(:,1)),t)
    results = [];
    results_time = [];
    results_trivial = [];
    
    for i = 1:n_graphs
        disp("Left: " + (n_graphs-i) + ", Size: " + size(G.Nodes, 1))
        test2do = n_tests(size(G.Nodes, 1), fract_targ, fract_dist, cluster_size, cluster_size);
        disp("test2do: " + test2do)
        for j = 1:test2do
            [results(end+1), results_time(end+1), results_trivial(edn+1)] = decouple(G, fract_targ, fract_dist, hubs(t,:), hubs(d,:));
        end
    end
    graph_name = t + "-" + d;% convertCharsToStrings(tags{tag});
    subplot(1,3,1);
    hold on
    add2boxchart(results, graph_name, "Cost", "Cost [-]", "Hubs")
    hold off
    subplot(1,3,2);
    hold on
    add2boxchart(results_time, graph_name, "Electrical network" + newline + "Runtime", "Time [s]", "Hubs")
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
saveas(gcf, "figures_new/Electrical network hubs.fig")

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