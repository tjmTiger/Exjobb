function  G = sfg(n, alpha, beta, gamma, delta_in, delta_out, seed)
% SFG
%    Description:
%       Returns a scale-free directed graph.
%       This implementation is based on a python implementaion from 
%       https://networkx.github.io/documentation/latest/reference/generators.html
%       Note: Sum of "alpha", "beta", and "gamma" must be 1.
%    Output Arguments:
%       G      : generated random digraph
%    Input Arguments:
%       n      : integer, number of nodes in graph
%       alpha  : float, probability for adding a new node connected to an
%                existing node chosen randomly according to the in-degree distribution.
%       beta   : float, probability for adding an edge between two existing nodes.
%                One existing node is chosen randomly according the in-degree 
%                distribution and the other chosen randomly according to the
%                out-degree distribution.     
%       gamma  : float, probability for adding a new node connected to an existing node
%                chosen randomly according to the out-degree distribution.
%       delta_in : float, bias for choosing ndoes from in-degree distribution.
%       delta_out: float, bias for choosing ndoes from out-degree distribution.
%       seed   : integer, seed for rng.
%    References:
%           [1] B. Bollob?s, C. Borgs, J. Chayes, and O. Riordan,
%               Directed scale-free graphs,
%               Proceedings of the fourteenth annual ACM-SIAM Symposium on
%               Discrete Algorithms, 132--139, 2003.
arguments
    n {mustBeNumeric}
    alpha single
    beta single
    gamma single
    delta_in {mustBeNumeric}
    delta_out {mustBeNumeric}
    seed {mustBeNumeric}
end
rng(seed)

% Chack if parameters are ok
if alpha < 0
    disp('alpha must be >= 0.')
end
if beta < 0
    disp('beta must be >= 0.')
end
if gamma < 0
    disp('gamma must be >= 0.')
end

if abs(alpha+beta+gamma - 1)>1e-10
    beta = beta - (alpha+beta+gamma - 1);
    disp("alpha+beta+gamma must equal 1. Beta is set to " + beta + " to compensate.")
end

% Initial graph (ring lattice)
G = speye(9);
G = [[sparse(1,9) 1];[G sparse(9,1)]];

% Add nodes untill graph has size n.
while size(G,1) < n
    r = rand();
    n_now = size(G,1);
    % random choice in alpha, beta, gamma ranges
    % alpha
    if r<alpha
        % add new node v
        G = [G sparse(n_now,1);sparse(1,n_now+1)];
        v = n_now+1;
        % choose w according to in-degree and delta_in
        w = choose_node(G, sum(G,1),delta_in);
    % beta
    elseif r < alpha+beta
        % choose v according to out-degree and delta_out
        v = choose_node(G, sum(G,2),delta_out);
        % choose w according to in-degree and delta_in
        w = choose_node(G, sum(G,1),delta_in);
    % gamma
    else
        % choose v according to out-degree and delta_out
        v = choose_node(G, sum(G,2),delta_out);
        % add new node w
        G = [G sparse(n_now,1);sparse(1,n_now+1)];
        w = n_now+1;
    end
    G(v,w) = 1; % add previously defined edge

end
% remove self-loops (not sure if there can even be self loops to begin with)
G = G - diag(diag(G));

G = digraph(G);
end

function i = choose_node(G,distribution,delta)
    % no idea how this part works
    cumsum_ = cumsum(distribution+delta);
    cumsum_ = cumsum_./cumsum_(end);
    r=rand();
    i = find(r<cumsum_,1);
end