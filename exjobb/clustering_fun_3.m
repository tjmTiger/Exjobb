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

% --- all-pairs shortest paths ---
D = distances(G);      % n x n matrix, Inf if unreachable

%% ----------------- Constraint building (vectorized) -----------------

rows = []; cols = []; vals = [];
bineq = [];

rowIdx = 1;

% 1) disjointness: x_i + y_i <= 1
rows = [rows, 1:n, 1:n];
cols = [1:n, n+1:2*n];
vals = [ones(1,n), ones(1,n)];
bineq = [bineq; ones(n,1)];
rowIdx = n+1;

% 2) intra-cluster (strong ε): forbid if not mutually ≤ ε
feasibleIntra = (D <= epsilon) & (D' <= epsilon); % symmetric
mask = triu(~feasibleIntra,1);   % only upper triangle, i<j
[i,j] = find(mask);

% forbid both in A
rows = [rows, rowIdx:rowIdx+numel(i)-1, rowIdx:rowIdx+numel(i)-1];
cols = [cols, i', j'];
vals = [vals, ones(1,numel(i)), ones(1,numel(i))];
bineq = [bineq; ones(numel(i),1)];
rowIdx = rowIdx + numel(i);

% forbid both in B
rows = [rows, rowIdx:rowIdx+numel(i)-1, rowIdx:rowIdx+numel(i)-1];
cols = [cols, (n+i)', (n+j)'];
vals = [vals, ones(1,numel(i)), ones(1,numel(i))];
bineq = [bineq; ones(numel(i),1)];
rowIdx = rowIdx + numel(i);

% 3) asymmetric inter-cluster separation: require d(i->j) >= gamma
mask = (D < gamma) & ~eye(n); % disallow A_i with B_j if dist < gamma
[i,j] = find(mask);

rows = [rows, rowIdx:rowIdx+numel(i)-1, rowIdx:rowIdx+numel(i)-1];
cols = [cols, i', (n+j)'];
vals = [vals, ones(1,numel(i)), ones(1,numel(i))];
bineq = [bineq; ones(numel(i),1)];
rowIdx = rowIdx + numel(i);

% build sparse Aineq
nConstr = rowIdx-1;
Aineq = sparse(rows, cols, vals, nConstr, 2*n);

%% ----------------- Solve ILP -----------------
Nvars = 2*n;
intcon = 1:Nvars;
lb = zeros(Nvars,1);
ub = ones(Nvars,1);

f = -ones(Nvars,1);   % maximize |A|+|B|

opts = optimoptions('intlinprog','Display','off'); % strict solver
% opts = optimoptions('intlinprog', ... % more relaxed solver
%     'Display','none', ...
%     'CutGeneration','none', ...         % fewer cuts for speed
%     'Heuristics','advanced', ...
%     'RelativeGapTolerance',0.05, ...    % allow 5% optimality gap
%     'MaxTime',30);                      % 30s limit

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
