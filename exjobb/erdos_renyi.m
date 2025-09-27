function [Gg,n,m] = erdos_renyi(n,p,seed)
% ERDOS_RENYI
%    Description:
%        Function generates a directed Erdos-Renyi random graph without self-loops.
%        Note: This code only generate approximately Erdos-Renyi random graph. 
%        Since Erdos-Renyi model only consider the undirected, non-self-loop
%        graphs. However, this code creates a directed graph, deletes
%        self-loops and only returns largest connected component.
%        
%        However, when the graph size n and probability p are large enough, the generated graph would
%        be almost same as Erdos-Renyi Model with directed edges.
%    Output Arguments:
%        Gg : generated random digraph
%        n  : graph size, number of nodes, |V|
%        m  : graph size, number of edges, |E|
%    Input Arguments:
%        n : integer, graph size, number of nodes, |V|
%        p : float, the probability p in the second definition of Erdos-Renyi
%            model. Note: this algorithm works only well for p << 1! (see sprand())
%        seed: seed for rng.
arguments
    n {mustBeNumeric}
    p single
    seed {mustBeNumeric}
end
rng(seed);
G = spones(sprand(n, n, p));  % sparse matrix with edges probability p, spones makes it unweighted
G = G - diag(diag(G));        % remove self-loops (set diagonal to zeros)
Ag = full(G);

% Get largest connected element
[bin, binsize] = conncomp(graph(sparse(Ag + Ag')));
n_comp = length(binsize);
ind_comp = bin;

if n_comp > 1 % If disconnected, make array of disconnected parts of our graph.
    for i = 1:n_comp
        A_dir_com{i} = Ag(ind_comp==i,ind_comp==i);
    end
else % if connected, put just out graph in the array.
    A_dir_com{1} = Ag;
end

% sort disconnected parts by size.
[~,ind_sort_com] = sort(binsize,'descend'); 

% Get largest graph according to ind_sort_com created earlier
Agg = A_dir_com{ind_sort_com(1)};
n = size(Agg,1);

% if(numel(A_dir_com)>1)
%     disp("WARNING: graph disconnected, largest component was used instead, with size: " + n)
% end

if nargout>2
    m = nnz(G); % Note: number of non zero elements
end
Gg = digraph(Agg');
end