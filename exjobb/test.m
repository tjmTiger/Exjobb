clear all;
clear;
clc;
G = ER_Graph(10, 0.5, 0);
% G = WattsStrogatz(10, 2, 0);
% G = SFG_dir(100, 0.3, 0.5, 0.2, 0);

figure(1);
%plot(G)
plot(G,'NodeColor','k','Layout','circle');

% load 'konect.mat';
% for i = 1:10
%     G = format_konect(data{i});
%     figure(i);
%     plot(G)
% end
% 
% load 'netzschleuder.mat';
% for j = 1:10
%     G = format_netz(data{j});
%     figure(i+j);
%     plot(G)
% end


%-------------%
%  Functions  %
%-------------%

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