function [Gg,n,m] = ER_Graph(n,p,s,seed,format,opt)
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
%        format:
%        opt:
switch nargin
    case 2
        seed = 0;
        format = 1;
        verbose = false;
    case 3
        format = 1;
        verbose = false;
    case 4
        verbose = false;
    otherwise
        disp('input argument invalid')
    
end
% rng(seed);
G = spones(triu(sprand(n,n,p),1)); % Note: random upper triangle of a binary (sparse) matrix (with density p). middle and lower triangles are zeros
if nargout>2
    m = nnz(G); % Note: number of non zero elements
end
if format==1
    G = G + G'; % Note: symetrical matrix?, this format shows entire matrix (not just half)?    
end

Ag = full(G); % Note: converts to normal matrix

for i = 1:size(Ag,1)
    for j = i+1 : size(Ag,1)
        if Ag(i,j) ~= 0
            if rand >= 0.5
                Ag(i,j) = 0;
            else
                Ag(j,i) = 0;
            end
        end
    end
end

A_undir = Ag+Ag';

[bin, binsize] = conncomp(graph(sparse(abs(A_undir))));
n_comp = length(binsize);
ind_comp = bin;

if n_comp > 1
%     s = sprintf('NB: undirected graph has %g connected components \n',n_comp);
%     disp(s);
    for i = 1:n_comp
        A_undir_com{i} = A_undir(ind_comp==i,ind_comp==i);
        % get the corresponding comp in the direct graph
        A_dir_com{i} = Ag(ind_comp==i,ind_comp==i);
%         s = sprintf('%g-th subgraph:   %g edges,  %g nodes \n',i, nnz(A_undir_com{i})/2, size(A_undir_com{i},1));
%         disp(s);
    end
else
    A_undir_com{1} = A_undir;
    A_dir_com{1} = Ag;
end

for i = 1:numel(A_undir_com)
    size_com(i) = size(A_undir_com{i},1);
end

[~,ind_sort_com] = sort(size_com,'descend');

Agg = A_dir_com{ind_sort_com(1)};
n = size(Agg,1);

Gg = digraph(Agg');

end