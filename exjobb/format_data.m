clear;
close all;
clc;

formated_data = dictionary("start", {1});
load 'netzschleuder.mat';
for j = 1:numel(info)
    tag = convertCharsToStrings(info{j}.tag);
    G = format_netz(data{j});
    if isKey(formated_data, tag)
        tmp = formated_data(tag);
        tmp = tmp{1,:};
        formated_data(tag) = {{tmp{1,:} {G}}}; % {"decouple(G, fract_targ, fract_dist)"};
    else
        formated_data(tag) = {{{G}}}; % {"decouple(G, fract_targ, fract_dist)"};
    end
    % results(1, end+1) = decouple(G, fract_targ, fract_dist);
end
formated_data("start") = [];

load 'konect.mat';
for j = 1:numel(info)
    tag  = convertCharsToStrings(info{j}.tag);
    if tag ~= ""
        G = format_konect(data{j});
        if isKey(formated_data, tag)
            tmp = formated_data(tag);
            tmp = tmp{1,:};
            formated_data(tag) = {{tmp{1,:} {G}}}; % {"decouple(G, fract_targ, fract_dist)"};
        else
            formated_data(tag) = {{{G}}}; % {"decouple(G, fract_targ, fract_dist)"};
        end
    end
end
save("formated_data.mat", "formated_data");













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