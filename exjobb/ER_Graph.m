function [Gg,n,m] = ER_Graph(n,p,s,seed) % cleaned up version
%    Description:
%        this function create Erdos-Renyi random Graph*
%        *This code only generate approximately Erdos-Renyi Random Graph. 
%        Since Erdos-Renyi Model only consider the undirected, non-self-loop
%        graphs. However, this code would firstly create a directed graph with,
%        self-loops. And then transform the directed graph into undirected simply
%        by ignore the upper triangular adjacency matrix and delete the self-loops
%        
%        However, when the graph size n is large enough, the generated graph would
%        be approximately similar to the expected Erdos-Renyi Model.
%    Output Arguments:
%        G : generated random graph
%        n : graph size, number of vertexes, |V|
%        m : graph size, number of edges, |E|
%    Input Arguments:
%        n : graph size, number of vertexes, |V|
%        p : the probability p of the second definition of Erdos-Renyi model.
%        s : 1 or 0 (if connceted component graph has nn=n or nn <= n)
%        seed: seed of the function. 

switch nargin
    case 2
        format = 1;
        verbose = false;
    case 3
        seed = 0;
        rng(seed);
        format = 1;
        verbose = false;
    case 4
        rng(seed)
        format = 1;
        verbose = false;
    otherwise
        disp('input argument invalid')
    
end
G = spones(triu(sprand(n,n,p),1)); % Note: random upper triangle of a binary (sparse) matrix (with density p). middle and lower triangles are zeros
if nargout>2
    m = nnz(G); % Note: number of non zero elements
end
if format==1
    G = G + G'; % Note: symetrical matrix  
end

Ag = full(G); % Note: converts to normal matrix

[bin, binsize] = conncomp(graph(sparse(abs(Ag))));
n_comp = length(binsize);
ind_comp = bin;

if n_comp > 1
    for i = 1:n_comp % Make array of graphs in out graph
        A_dir_com{i} = Ag(ind_comp==i,ind_comp==i);
    end
else
    A_dir_com{1} = Ag;
end

[~,ind_sort_com] = sort(binsize,'descend'); 


Agg = A_dir_com{ind_sort_com(1)};% Get largest graph according to ind_sort_com created earlier
n = size(Agg,1);

if(numel(A_dir_com)>1)
    disp("WARNING: graph disconnected, largest component was used instead, with size: " + n)
end

for i = 1:size(Agg,1) % Makes Ag directional
    for j = i+1 : size(Agg,1)
        if Agg(i,j) ~= 0
            if rand >= 0.5
                Agg(i,j) = 0;
            else
                Agg(j,i) = 0;
            end
        end
    end
end

Gg = digraph(Agg');
end