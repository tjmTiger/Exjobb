function h = WattsStrogatz(N,K,beta, seed)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

switch nargin
    case 3
        disp("no seed")
    case 4
        rng(seed);
    otherwise
        disp('input argument invalid')
end

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K); % Note: array with once in first row, copunting up on every next row
t = s + repmat(1:K,N,1); % Note: s made into diagonal matrix instead, 2 in left upper corner, counting up diagonally
t = mod(t-1,N)+1; % Note: mod returns remainder after division of t-1 by N => makes N to largest number in the matrix, after that reapeating from 1

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = digraph(s,t);

% h = graph(s,t);

% Aa = full(adjacency(h));
% 
% for i = 1:size(Aa,1) % make graph directional
%     for j = i+1 : size(Aa,1)
%         if Aa(i,j) ~= 0
%             if rand >= 0.5
%                 Aa(i,j) = 0;
%             else
%                 Aa(j,i) = 0;
%             end
%         end
%     end
% end
% h = digraph(Aa');


end