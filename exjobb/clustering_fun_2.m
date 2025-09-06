%% ILP for two disjoint epsilon-clubs with strong separation (MATLAB)
% Usage: run as-is. Produces a random test adjacency matrix A (20x20),
% builds digraph G = digraph(A'), computes all-pairs shortest paths,
% and solves ILP to maximize |A|+|B| under the constraints described.

clear; clc;

t_start = tic;
%% ILP for two disjoint epsilon-clubs with asymmetric separation (MATLAB)
epsilon = 4;           % intra-cluster diameter bound
gamma = 2;             % inter-cluster separation bound
rng(1);                % reproducible

% test graph
n = 20;                % number of nodes
% --- test adjacency matrix A and digraph G = digraph(A') ---
p = 0.1;              % edge probability for random digraph
A = rand(n) < p;       % random adjacency
A(1:n+1:end) = 0;      % no self-loops
G = digraph(A');       % note: digraph(A')

% % get network from data file
% tag = 12; % 11 = Technological network; 12 = Electrical network;
% load formated_data.mat;
% tags = keys(formated_data);
% val = values(formated_data);
% disp(tags{tag})
% G = val{tag}{1}{1};
% n = numnodes(G);


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

% --- all-pairs shortest paths ---
D = distances(G);      % n x n matrix, Inf if unreachable

%% ----------------- Build constraints efficiently -----------------
rows = []; cols = []; vals = [];
bineq = [];

rowIdx = 1;

% 1) disjointness: x_i + y_i <= 1
for i=1:n
    rows = [rows rowIdx rowIdx];
    cols = [cols i n+i];
    vals = [vals 1 1];
    bineq(rowIdx,1) = 1;
    rowIdx = rowIdx + 1;
end

% 2) intra-cluster (strong ε): forbid if not mutually ≤ ε
for i=1:n-1
    for j=i+1:n
        dij = D(i,j); dji = D(j,i);
        if ~((dij <= epsilon) && (dji <= epsilon))
            % forbid both in A
            rows = [rows rowIdx rowIdx];
            cols = [cols i j];
            vals = [vals 1 1];
            bineq(rowIdx,1) = 1;
            rowIdx = rowIdx + 1;

            % forbid both in B
            rows = [rows rowIdx rowIdx];
            cols = [cols n+i n+j];
            vals = [vals 1 1];
            bineq(rowIdx,1) = 1;
            rowIdx = rowIdx + 1;
        end
    end
end

% 3) asymmetric inter-cluster separation: require d(i->j) >= gamma
for i=1:n
    for j=1:n
        if i~=j && D(i,j) < gamma
            rows = [rows rowIdx rowIdx];
            cols = [cols i n+j];
            vals = [vals 1 1];
            bineq(rowIdx,1) = 1;
            rowIdx = rowIdx + 1;
        end
    end
end

% build sparse Aineq
nConstr = rowIdx-1;
Aineq = sparse(rows, cols, vals, nConstr, 2*n);

%% ----------------- Solve ILP -----------------
Nvars = 2*n;
intcon = 1:Nvars;
lb = zeros(Nvars,1);
ub = ones(Nvars,1);

f = -ones(Nvars,1);   % maximize |A|+|B|

opts = optimoptions('intlinprog','Display','off');
[z,fval,exitflag,output] = intlinprog(f,intcon,Aineq,bineq,[],[],lb,ub,opts);

x = round(z(1:n));
y = round(z(n+1:2*n));

Aset = find(x==1);
Bset = find(y==1);

%% ----------------- Report results -----------------
fprintf('n = %d, epsilon = %d, gamma = %d\n', n, epsilon, gamma);
fprintf('Solution: |A|=%d, |B|=%d, total=%d\n', numel(Aset), numel(Bset), numel(Aset)+numel(Bset));
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

disp(toc(t_start));

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
