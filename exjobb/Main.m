clear;
close all;
clc;

n_graphs = 100; % number of graphs
% fract_targ = 0.1;
% fract_dist = 0.1;

%-----------------------------------------------%
%                                               %
%                  Erdos Renyi                  %
%                                               %
%-----------------------------------------------%

% Distrubance and target node fractions

results_all = [];
results_time_all = [];
fract_T_D = 0.05:0.05:0.95;
for j = fract_T_D
    fract_targ = j;
    fract_dist = j;
    results = zeros(1, n_graphs);
    results_time = zeros(1, n_graphs);
    for i = 1:n_graphs
        n = 100;
        p = 2*(log10(n)/n);
        G = ER_Graph(n, p, 1, i);
        t_start = tic;
        results(i) = decouple(G, fract_targ, fract_dist);
        results_time(i) = toc(t_start);
    end
    results_all(end+1,:) = results;
    results_time_all(end+1,:) = results_time;
end

figure();
hold on;

for r = 1:size(results_all, 1)
    add2boxchart(results_all(r,:), "frac: " + fract_T_D(r), "Erdos Renyi, T & D fractions")
end

hold off;

figure();
hold on;

for r = 1:size(results_time_all, 1)
    add2boxchart(results_time_all(r,:), "frac: " + fract_T_D(r), "Erdos Renyi, T & D fractions, elapsed time", "time [s]")
end

hold off;

% Size

% Edge probability

%-----------------------------------------------%
%                                               %
%                 Watts Strogatz                %
%                                               %
%-----------------------------------------------%
%
% results = zeros(1, n_graphs);
% for i = 1:n_graphs
%     G = WattsStrogatz(50, 2, 0.2);
%     results(i) = decouple(G, fract_targ, fract_dist);
% end
% add2boxchart(results, "Watts Strogatz")
% 
%-----------------------------------------------%
%                                               %
%                   Scale Free                  %
%                                               %
%-----------------------------------------------%
%
% results = zeros(1, n_graphs);
% for i = 1:n_graphs
%     G = SFG_dir(50, 0.2, 0.2, 0.6);
%     results(i) = decouple(G, fract_targ, fract_dist);
% end
% add2boxchart(results, "SFG")
% 
%-----------------------------------------------%
%                                               %
%                  Real Networks                %
%                                               %
%-----------------------------------------------%
%
% results = [];
% load 'netzschleuder.mat';
% tag = "None";
% for j = 1:numel(info)
%     new_tag = convertCharsToStrings(info{j}.tag);
%     if ~strcmp(new_tag, tag) && ~strcmp(tag, "None")
%         add2boxchart(results, tag);
%         results = [];
%     end
%     tag = new_tag;
%     G = format_netz(data{j});
%     disp(size(G.Nodes));
%     results(1, end+1) = decouple(G, fract_targ, fract_dist);
% end
% add2boxchart(results, tag);
%
% results = [];
% load 'konect.mat';
% tag = "None";s
% j = 1;
% while numel(info) >= j
%     tag  = convertCharsToStrings(info{j}.tag);
%     disp(tag)
%     i = 1;
%     while numel(info) >= i
%         if tag == convertCharsToStrings(info{i}.tag)
%             G = format_konect(data{i});
%             % G = rmedge(G, 1:numnodes(G), 1:numnodes(G)); %     remove self loops
%             disp(size(G.Nodes));
%             results(1, end+1) = decouple(G, fract_targ, fract_dist);
%             data(i) = [];
%             info(i) = [];
%         else
%             i=i+1;
%         end
%     end
%     add2boxchart(results, tag);
%     results = [];
%     tag = "None";
%     j = j+1;
% end

%-----------------------------------------------%
%                                               %
%                   Functions                   %
%                                               %
%-----------------------------------------------%

function add2boxchart(results, test_name, title_name, ylabel_name)
    boxchart(categorical(1:numel(results), 1:numel(results), repmat(test_name, 1, numel(results))), results)
    if nargin >= 3
        title(title_name)
    end
    if nargin >= 4
        ylabel(ylabel_name)
    else
        ylabel('Cost')
    end
    xlabel('Tests')
end

% function add2bar(results)
%     figure();
%     hold on;
%     elements = unique(results);
%     for x = elements
%         y = sum(ismember(results, x));
%         bar(x, y, 'FaceColor',[0.3 0.6 0.6],'EdgeColor',[0 0 0],'LineWidth',0.5)
%     end
%     hold off;
%     xlabel("Cost");
%     ylabel("Num. of graphs")
% end

function G = format_konect(G)
    if ~is_directed(G)
        G = G+G';
    end
    G = get_largest(G);
    G = digraph(G');
end

function G = format_netz(G)
    G = full(G);
    if ~is_directed(G)
        G = G+G';
    end
    G = get_largest(G);
    G = digraph(G');
end

function dir = is_directed(G)
    lower = adjacency(graph(G, 'lower'));
    upper = adjacency(graph(G, 'upper'));
    dir = isempty(find(lower, 1)) == isempty(find(upper, 1));
end

function G = get_largest(G)
    [bin, binsize] = conncomp(graph(sparse(abs(G+G')), 'upper'));
    n_comp = length(binsize);
    ind_comp = bin;
    
    if n_comp > 1
        for i = 1:n_comp % Make array of graphs in out graph
            A_dir_com{i} = G(ind_comp==i,ind_comp==i);
        end
    else
        A_dir_com{1} = G;
    end
    
    [~,ind_sort_com] = sort(binsize,'descend'); 
    
    G = A_dir_com{ind_sort_com(1)};% Get largest graph according to ind_sort_com created earlier
end