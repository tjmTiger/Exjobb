function Gg = watts_strogatz(n,k_half,p_ws, seed)
% WATTS_STROGATZ
%    Description:
%       Returns a Watts-Strogatz model graph with n nodes, n*k edges,
%       mean node degree k = k_half*2, and rewiring probability p_ws.
%       p_ws ∈ [0, 1] where p_ws = 0 is a ring lattice, and p_ws = 1 is a random graph.
%    Output Arguments:
%       Gg      : generated random digraph
%    Input Arguments:
%       n       :desired graph size (number of nodes)
%       k       : integer, avergae degree k.
%       p_ws    : float, probability for rewiring.
%       seed    : integer, seed for rng.
arguments
    n {mustBeNumeric}
    k_half {mustBeNumeric}
    p_ws single
    seed {mustBeNumeric}
end
% Connect each node to its k neighbors. This constructs indices for a ring lattice.
rng(seed);
s = repelem((1:n)',1,k_half); % note: array with once in first row, copunting up on every next row
t = s + repmat(1:k_half,n,1); % note: s made into diagonal matrix instead, 2 in left upper corner, counting up diagonally
t = mod(t-1,n)+1; % note: mod returns remainder after division of t-1 by n => makes n to largest number in the matrix, after that repeating from 1

% Rewire the target node of each edge with probability p_ws
for source=1:n    
    switchEdge = rand(k_half, 1) < p_ws; % defines if edge will be switched
    
    newTargets = rand(n, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

Gg = digraph(s,t);
end