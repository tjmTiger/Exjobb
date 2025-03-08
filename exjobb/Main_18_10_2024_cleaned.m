% cleaned up
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

G = ER_Graph(200, 0.5, 0);
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
A = A.*randn(N,N);
G = digraph(A');

[V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair_backprop_new(G, D, T); % DDDF


k = 1; % randi([1,numel(V_out)]); % select a (C,A,B)-pair
COM = compensator_ddp(G, D, T, V_out{k}, V_in{k}, S_Min{k}, Z_Max{k})

syms t real
t = 1;
expo = round(COM.T_e*expm(sparse(COM.A_e)*t)*COM.D_e) % ???

V_out = submincutDDSF_final2(G,D,T,'V_out');
[V_out_cell, G_cell] = mincutDDSF_all(G,D,T,V_out,'V_out','all'); % find optimal V_in & V_out
[V_out_cell2, G_cell2] = mincutDDSF_all(G,D,T,V_out_cell,'V_out','all',G_cell{end});
[V_out_cell3, G_cell3] = mincutDDSF_all(G,D,T,[V_out_cell, V_out_cell2],'V_out','all',G_cell2{end});
[V_out_cell4, G_cell4] = mincutDDSF_all(G,D,T,[V_out_cell, V_out_cell2, V_out_cell3],'V_out','all',G_cell3{end});
[V_out_cell5, G_cell5] = mincutDDSF_all(G,D,T,[V_out_cell, V_out_cell2, V_out_cell3, V_out_cell4],'V_out','all',G_cell4{end});

GG = minreal(tf(ss(COM.A_e,COM.D_e,COM.T_e,zeros(length(T),length(D)))));

k = 1;
COM = compensator_ddp(G, D, T, V_out, V_in{k}, S_Min{k}, Z_Max{k}) % typ som framkoppling?

syms t real
t = 1;
expo = COM.T_e*expm(COM.A_e*t)*COM.D_e
% DDSF output nodes -> V_out; DDOF output nodes -> Vout
V_in = submincutDDSF_final2(G,D,T,'V_in');
V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');
[Vin_opt, Vout, C1] = constrained_optimal_solution(G,D,T,V_in_all,'V_out'); % DDOF, Given the list of optima V_in (or V_out), computes all V_out-V_in
V_out = submincutDDSF_final2(G,D,T,'V_out');
V_out_all = mincutDDSF_all(G,D,T,V_out,'V_out','all');
[Vin, Vout_opt, C2] = constrained_optimal_solution(G,D,T,V_out_all,'V_in'); % DDOF
[V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(Vin_opt, Vout, C1, Vin, Vout_opt, C2); % DDOF, Given the list of optimal lists of pairs V_in-V_out provided by constrained_optimal_solution.m find min C 

if ~isempty(V_in_best)

    for x = 1:numel(V_in_best)
        figure;
        p = plot(G,'b');
        title('$\mathcal{G}$')
        nodeColors = 1 * ones(N, 1); % Default to value 3 (Yellow)
        nodeColors(V_in_best{x}) = 3;
        nodeColors(V_out_best{x}) = 2;
        p.NodeCData = nodeColors;
        colormap(jet); % Use the 'jet' colormap
        p.MarkerSize = 8; % Increase or decrease the size of the nodes
        hold on; % Hold on to the current plot
        legendEntries = {'Input nodes','Output nodes','Other nodes'};
        hRed = scatter(nan, nan, 100, 'r', 'filled'); % Placeholder for red nodes
        hGreen = scatter(nan, nan, 100, 'g', 'filled'); % Placeholder for green nodes
        hBlue = scatter(nan, nan, 100, 'b', 'filled'); % Placeholder for blue nodes
        legend([hRed, hGreen, hBlue], legendEntries, 'Location', 'best');
        hold off; % Release the hold on the current plot
    end

end
