function [cost, results_time, trivial_solutions] = decouple(G, fract_targ, fract_dist)
% Note: frac_targ + frac_dist <= 1

if fract_targ + fract_dist > 1
    error('input argument invalid, fract_targ + fract_dist must be lesst than 1')
end

% cleaned up
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

N = numnodes(G);

% fract_targ = 0.001;
% fract_dist = 0.001;
T = [];

while isempty(T) % add targets and disturbances
    n_dist = ceil(fract_dist*N);
    n_targ = ceil(fract_targ*N);
    D = sort(randsample(N, n_dist));
    T = sort(randsample(setdiff(1:N', D), n_targ))';
    % T = sort(setdiff(randsample(N, n_targ),D));
end

n_dist = length(D);
n_targ = length(T);

V_in_initial = []; % control on targets if those are directly connected to a disturbance
population = setdiff(setdiff(1:N, T), D);

for i = 1:n_targ % decouple targets from disturbances and other targets
    a = 0;
    [sout, ~] = findedge(G,inedges(G,T(i))); % find nodes with outgoing edges onto targets§
    if ~isempty(intersect(sout,D)) % check if those are disturbance nodes
        DT = intersect(sout,D);
        L1 = length(DT);
        V_in_initial = [V_in_initial, T(i)]; % if yes put the control on the target nodes
        G = rmedge(G,DT,T(i)); % remove incoming edges of targets
    end
    if a
        [sout_2, ~] = findedge(G,inedges(G,T(i)));
        if isempty(sout_2)
            s = randsample(population, min(length(population),ceil(L1)));
            G = addedge(G, s, T(i), 1);
        end
    end
end

A = full(adjacency(G))';
% A = A.*randn(N,N); % random weights
G = digraph(A');

% V_out = submincutDDSF_final2(G,D,T,'V_out');
% V_out_all = mincutDDSF_all(G,D,T,V_out,'V_out','all');
t_start = tic;
V_in = submincutDDSF_final2(G,D,T,'V_in');
results_time = toc(t_start);
V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');
% disp("n_dist = " + n_dist)
% disp("n_targ = " + n_targ)
% disp(V_in)
% disp(V_out)
% disp(D)
% disp(T)

%check and display how many control nodes are placed on target nodes
v_in_on_T = numel(V_in);
for v_in = V_in_all
    v_in = v_in{1};
    [C,~,ic] = unique([v_in T']);
    a_counts = accumarray(ic,1);
    v_in_on_T_next = sum(a_counts(:,1)~=1);
    if v_in_on_T > v_in_on_T_next
            v_in_on_T = v_in_on_T_next;
    end
end

trivial_solutions = v_in_on_T/numel(V_in);
% if v_in_on_T ~= 0
%     disp("WARNING: There are " + v_in_on_T + " (out of " + numel(V_in) + ") measured nodes placed on target nodes.")
% end

% [C,~,ic] = unique([V_in T']);
% a_counts = accumarray(ic,1);
% 
% if ~isequal(a_counts, ones(size(a_counts)))
%     disp("WARNING: There are " + sum(a_counts(:,1)~=1) + " (out of " + numel(V_in) + ") measured nodes placed on target nodes.")
% end

% disp("numebr of V_in = " + numel(V_in))
% disp("number of V_out = " + numel(V_out))
% disp("number of V_in + V_out = " + (numel(V_in) + numel(V_out)))
% disp("-----------------")
cost = ( numel(V_in)) / ( n_targ + n_dist );
% cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );