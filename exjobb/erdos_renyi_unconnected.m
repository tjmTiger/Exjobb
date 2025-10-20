function [Gg,n,m] = erdos_renyi_unconnected(n,p,seed)
% ERDOS_RENYI
%    Description:
%        Function generates a directed Erdos-Renyi random graph without self-loops.
%        Note: This code only generate approximately Erdos-Renyi random graph. 
%        Since Erdos-Renyi model only consider the undirected, non-self-loop
%        graphs. However, this code creates a directed graph and deletes
%        self-loops.
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

if nargout>2
    m = nnz(G); % Note: number of non zero elements
end
Gg = digraph(Ag');
end