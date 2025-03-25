clear;
close all;
clc;
% cleaned up
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

G = ER_Graph(10, 0.6, 0);
plot(G)
N = numnodes(G);

fract_targ = 0.001;
fract_dist = 0.001;
T = [];

while isempty(T) % add targets and disturbances?
    n_dist = ceil(fract_dist*N);
    n_targ = ceil(fract_targ*N);
    D = sort(randsample(N, n_dist));
    T = sort(setdiff(randsample(N, n_targ),D));
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


figure(); % After edge removal by action of V_in_initial on targets directly connected to disturbances
p = plot(G,'b');
title('$\mathcal{G}$')
nodeColors = 1 * ones(N, 1); % Default to value 3 (Yellow)
nodeColors(T) = 2;
nodeColors(D) = 3;
p.NodeCData = nodeColors;
colormap(jet); % Use the 'jet' colormap
p.MarkerSize = 8; % Increase or decrease the size of the nodes
hold on; % Hold on to the current plot
legendEntries = {'Disturbance', 'Target', 'Other nodes'};
hRed = scatter(nan, nan, 100, 'r', 'filled'); % Placeholder for red nodes
hGreen = scatter(nan, nan, 100, 'g', 'filled'); % Placeholder for green nodes
hYellow = scatter(nan, nan, 100, 'b', 'filled'); % Placeholder for yellow nodes
legend([hRed, hGreen, hYellow], legendEntries, 'Location', 'best');
hold off; % Release the hold on the current plot



A = full(adjacency(G))';
% A = A.*randn(N,N); % random weights
G = digraph(A');

V_out = submincutDDSF_final2(G,D,T,'V_out');
V_out_all = mincutDDSF_all(G,D,T,V_out,'V_out','all');
V_in = submincutDDSF_final2(G,D,T,'V_in');
V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');

% disp(numel(V_in) + numel(V_out))
% disp("--------------")
% disp(V_in_all)
% disp(V_out_all)

min_cost = numel(V_in) + numel(V_out);

% T = table(V_in', V_out', cellfun(@numel, V_in)', cellfun(@numel, V_out_best)',  C.*ones(numel(S),1), S, 'VariableNames', {'V_in', 'V_out','# of Inputs','# of Outputs', 'Minimal Cost C', 'Constraint w.r.t.'});
% disp(T);