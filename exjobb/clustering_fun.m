%% ILP for two disjoint epsilon-clubs with strong separation (MATLAB)
% Usage: run as-is. Produces a random test adjacency matrix A (20x20),
% builds digraph G = digraph(A'), computes all-pairs shortest paths,
% and solves ILP to maximize |A|+|B| under the constraints described.

clear; clc;

%% ILP for two disjoint epsilon-clubs with asymmetric separation (MATLAB)

rng(1);                % reproducible
n = 20;                % number of nodes
epsilon = 4;           % intra-cluster diameter bound
gamma = 2;             % inter-cluster separation bound

% --- test adjacency matrix A and digraph G = digraph(A') ---
p = 0.1;              % edge probability for random digraph
A = rand(n) < p;       % random adjacency
A(1:n+1:end) = 0;      % no self-loops
G = digraph(A');       % note: digraph(A')
% adjacency matrix for digraph with 10 nodes
% n = 10;
% A = zeros(n);
% 
% edges = [1 3;
%          2 4;
%          2 3;
%          3 5;
%          4 6;
%          3 6;
%          5 7;
%          6 8;
%          6 7;
%          5 8;
%          7 9;
%          7,10;
%          8 10;
%          9 1];
% 
% for k = 1:size(edges,1)
%     i = edges(k,1);
%     j = edges(k,2);
%     A(i,j) = 1;
% end

% create digraph
% G = digraph(A);


% --- all-pairs shortest paths ---
D = distances(G);      % n x n matrix, Inf if unreachable

% --- ILP variables ---
% x(1:n) : node in A
% y(1:n) : node in B
Nvars = 2*n;
intcon = 1:Nvars;
lb = zeros(Nvars,1);
ub = ones(Nvars,1);

Aineq = [];
bineq = [];

% 1) disjointness: x_i + y_i <= 1
for i = 1:n
    row = zeros(1,Nvars);
    row(i) = 1; row(n+i) = 1;
    Aineq = [Aineq; row];
    bineq = [bineq; 1];
end

% 2) intra-cluster constraints (strong): forbid i,j both in A if either direction > epsilon
for i = 1:n-1
    for j = i+1:n
        dij = D(i,j); dji = D(j,i);
        if ~( (dij <= epsilon) && (dji <= epsilon) )
            % forbid both in A
            row = zeros(1,Nvars);
            row(i) = 1; row(j) = 1;
            Aineq = [Aineq; row];
            bineq = [bineq; 1];
            % forbid both in B
            row = zeros(1,Nvars);
            row(n+i) = 1; row(n+j) = 1;
            Aineq = [Aineq; row];
            bineq = [bineq; 1];
        end
    end
end

% 3) asymmetric inter-cluster separation:
%    only require d(i->j) >= gamma for all i in A, j in B
for i = 1:n
    for j = 1:n
        if i==j, continue; end
        if D(i,j) < gamma   % too close from i->j
            row = zeros(1,Nvars);
            row(i) = 1; row(n+j) = 1;
            Aineq = [Aineq; row];
            bineq = [bineq; 1];
        end
    end
end

% --- objective: maximize |A|+|B| ---
f = -ones(Nvars,1);

% --- solve ILP ---
opts = optimoptions('intlinprog','Display','off');
[z,fval,exitflag,output] = intlinprog(f,intcon,Aineq,bineq,[],[],lb,ub,opts);

x = round(z(1:n));
y = round(z(n+1:2*n));

Aset = find(x==1);
Bset = find(y==1);

fprintf('n = %d, epsilon = %d, gamma = %d\n', n, epsilon, gamma);
fprintf('Found solution: |A| = %d, |B| = %d, total = %d\n', numel(Aset), numel(Bset), numel(Aset)+numel(Bset));
fprintf('A nodes: '); fprintf('%d ', Aset); fprintf('\n');
fprintf('B nodes: '); fprintf('%d ', Bset); fprintf('\n');

% diagnostics
if ~isempty(Aset)
    maxDiamA = max(max(D(Aset,Aset)));
    fprintf('Max strong distance in A = %g\n', maxDiamA);
end
if ~isempty(Bset)
    maxDiamB = max(max(D(Bset,Bset)));
    fprintf('Max strong distance in B = %g\n', maxDiamB);
end
if ~isempty(Aset) && ~isempty(Bset)
    minAB = min(D(Aset,Bset),[],'all');
    fprintf('Min distance A->B = %g\n', minAB);
end

% visualization
figure;
h = plot(G,'Layout','force');
title('Digraph with A (green) and B (red)');
labelnode(h,1:n,1:n);
colors = repmat([0.8 0.8 0.8],n,1);
colors(Aset,:) = repmat([0 0.6 0],numel(Aset),1);
colors(Bset,:) = repmat([0.8 0 0],numel(Bset),1);
for k = 1:n
    highlight(h,k,'NodeColor',colors(k,:));
end
