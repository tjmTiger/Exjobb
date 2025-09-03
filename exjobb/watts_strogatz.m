function h = watts_strogatz(N,K_half,beta, seed)
% H = WattsStrogatz(N,K_half,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K_half edges, mean node degree 2*K_half, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

rng(seed)
% Connect each node to its K_half next and previous neighbors. This constructs
% indices for a ring lattice.
rng(seed);
s = repelem((1:N)',1,K_half); % Note: array with once in first row, copunting up on every next row
t = s + repmat(1:K_half,N,1); % Note: s made into diagonal matrix instead, 2 in left upper corner, counting up diagonally
t = mod(t-1,N)+1; % Note: mod returns remainder after division of t-1 by N => maK_halfes N to largest number in the matrix, after that reapeating from 1

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K_half, 1) < beta;
    
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
% for i = 1:size(Aa,1) % maK_halfe graph directional
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