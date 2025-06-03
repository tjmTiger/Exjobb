clear;
close all;
clc;

%-----------------------------------------------%
%                                               %
%                  Real Networks                %
%                                               %
%-----------------------------------------------%

fract_targ = 0.1;
fract_dist = 0.1;

load formated_data.mat;

tags = keys(formated_data);
val = values(formated_data);
figure();
results_all = [];
results_time_all = [];
results_trivial_all = [];
for tag = 1:length(tags)
    disp(tags{tag})

    n_graphs = length(val{tag});
    results = zeros(1, n_graphs);
    results_time = zeros(1, n_graphs);
    results_trivial = zeros(1, n_graphs);

    for i = 1:length(val{tag})
        G = val{tag}{i}{1};
        disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
        [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
    end
    graph_name = convertCharsToStrings(tags{tag});
    subplot(1,3,1);
    hold on
    add2boxchart(results, graph_name, "Cost", "Cost [-]", "Graph category")
    hold off
    subplot(1,3,2);
    hold on
    add2boxchart(results_time, graph_name, "Runtime", "Time [s]", "Graph category")
    hold off
    subplot(1,3,3);
    hold on
    add2boxchart(results_trivial, graph_name, "Trivial solutions", "Index [-]", "Graph category")
    hold off
end

load formated_data.mat;

tags = keys(formated_data);
val = values(formated_data);
figure();
for tag = 1:length(tags)
    disp(tags{tag})

    n_graphs = length(val{tag});
    results = zeros(1, n_graphs);
    results_time = zeros(1, n_graphs);
    results_trivial = zeros(1, n_graphs);

    for i = 1:length(val{tag})
        G = val{tag}{i}{1};
        disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
        results(i) = size(G.Nodes, 1);
    end
    graph_name = convertCharsToStrings(tags{tag});
    hold on
    add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
    hold off
end

load formated_data.mat;

tags = keys(formated_data);
val = values(formated_data);
% figure();
for tag = 1:length(tags)
    disp(tags{tag})

    n_graphs = length(val{tag});
    results = zeros(1, n_graphs);
    results_time = zeros(1, n_graphs);
    results_trivial = zeros(1, n_graphs);

    for i = 1:length(val{tag})
        G = val{tag}{i}{1};
        disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
    end
    % graph_name = convertCharsToStrings(tags{tag});
    % hold on
    % add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
    % hold off
end


load formated_data.mat;

tags = keys(formated_data);
val = values(formated_data);
for tag = 1:length(tags)
    disp(tags{tag})
    n_graphs = length(val{tag});

    graph_name = convertCharsToStrings(tags{tag});
    x = ceil(sqrt(n_graphs));
    figure();
    for i = 1:length(val{tag})
        G = val{tag}{i}{1};
        disp("Left: " + (length(val{tag})-i) + ", Size: " + size(G.Nodes, 1))
        subplot(x,x,i)
        plot(G)
    end
    sgtitle(graph_name)
    % hold on
    % add2boxchart(results, graph_name, "Real Networks, Sizes", "Size [psc]", "Graph category")
    % hold off
end

%%
%-----------------------------------------------%
%                                               %
%              Technological network            %
%                                               %
%-----------------------------------------------%
clear;
close all;
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


large_nodes = [];
for i = 1:length(hub_ranks)
    if hub_ranks(i) > max(hub_ranks)*0.7
        large_nodes(end+1,:) = [i,hub_ranks(i)];
    end
end
% disp(large_nodes)

color_n = 'r';
color_e = '#ff5555';
hub1 = [];
i_list = 482; % 650; % 482;
for j = 1:1
    for i = i_list
        p = predecessors(G,i);
        s = successors(G,i);
        highlight(fig,p,'NodeColor',color_n)
        highlight(fig,s,'NodeColor',color_n)
        highlight(fig,i,p,'EdgeColor',color_e)
        highlight(fig,i,s,'EdgeColor',color_e)
    end
    i_list = [p',s'];
    hub1 = [hub1, i_list];
end

color_n = 'b';
color_e = '#5577ff';
hub2 = [];
i_list = 445; % 292;
for j = 1:3
    for i = i_list
        p = predecessors(G,i);
        s = successors(G,i);
        highlight(fig,p,'NodeColor',color_n)
        highlight(fig,s,'NodeColor',color_n)
        highlight(fig,i,p,'EdgeColor',color_e)
        highlight(fig,i,s,'EdgeColor',color_e)
    end
    i_list = [p',s'];
    hub2 = [hub2, i_list];
end

color_n = 'g';
color_e = '#5aff55';
hub3 = [];
i_list = 347; % 99; % 445;
for j = 1:4
    for i = i_list
        p = predecessors(G,i);
        s = successors(G,i);
        highlight(fig,p,'NodeColor',color_n)
        highlight(fig,s,'NodeColor',color_n)
        highlight(fig,i,p,'EdgeColor',color_e)
        highlight(fig,i,s,'EdgeColor',color_e)
    end
    i_list = [p',s'];
    hub3 = [hub3, i_list];
end



figure();
fract_targ = 0.2;
fract_dist = 0.2;

n_graphs = 2;
results = zeros(1, n_graphs);
results_time = zeros(1, n_graphs);
results_trivial = zeros(1, n_graphs);

for i = 1:n_graphs
    disp("Left: " + (n_graphs-i) + ", Size: " + size(G.Nodes, 1))
    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist, hub1, hub2);
end
graph_name = "hub1-2"; % convertCharsToStrings(tags{tag});
subplot(1,3,1);
hold on
add2boxchart(results, graph_name, "Cost", "Cost [-]", "Graph category")
hold off
subplot(1,3,2);
hold on
add2boxchart(results_time, graph_name, "Runtime", "Time [s]", "Graph category")
hold off
subplot(1,3,3);
hold on
add2boxchart(results_trivial, graph_name, "Trivial solutions", "Index [-]", "Graph category")
hold off
ylim([0 1])

% ---------------------------------------- % 

results = zeros(1, n_graphs);
results_time = zeros(1, n_graphs);
results_trivial = zeros(1, n_graphs);

for i = 1:n_graphs
    disp("Left: " + (n_graphs-i) + ", Size: " + size(G.Nodes, 1))
    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist, hub2, hub3);
end
graph_name = "hub2-3";% convertCharsToStrings(tags{tag});
subplot(1,3,1);
hold on
add2boxchart(results, graph_name, "Cost", "Cost [-]", "Graph category")
hold off
subplot(1,3,2);
hold on
add2boxchart(results_time, graph_name, "Runtime", "Time [s]", "Graph category")
hold off
subplot(1,3,3);
hold on
add2boxchart(results_trivial, graph_name, "Trivial solutions", "Index [-]", "Graph category")
hold off
ylim([0 1])

% ---------------------------------------- % 

results = zeros(1, n_graphs);
results_time = zeros(1, n_graphs);
results_trivial = zeros(1, n_graphs);

for i = 1:n_graphs
    disp("Left: " + (n_graphs-i) + ", Size: " + size(G.Nodes, 1))
    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist, hub1, hub3);
end
graph_name = "hub1-3";% convertCharsToStrings(tags{tag});
subplot(1,3,1);
hold on
add2boxchart(results, graph_name, "Cost", "Cost [-]", "Graph category")
hold off
subplot(1,3,2);
hold on
add2boxchart(results_time, graph_name, "Runtime", "Time [s]", "Graph category")
hold off
subplot(1,3,3);
hold on
add2boxchart(results_trivial, graph_name, "Trivial solutions", "Index [-]", "Graph category")
hold off
ylim([0 1])

%%
%-----------------------------------------------%
%                                               %
%                Electrical network             %
%                                               %
%-----------------------------------------------%
% tag = 12

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