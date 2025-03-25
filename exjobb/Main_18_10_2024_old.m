clear;
close all;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

N = 200;
 [G,n] = ER_Graph(N,0.5,0);
 plot(G)
n = 105;
k = 205;
G = generate_directed_connected_graph_try(n,k);
A = full(adjacency(G))';

fract_targ = 0.001;
fract_dist = 0.001;
T = [];
N = size(A,1);

while isempty(T)
    n_dist = ceil(fract_dist*N); % randi([10, 50]);
    n_targ = ceil(fract_targ*N);
    D = sort(randsample(N, n_dist));
    T = sort(setdiff(randsample(N, n_targ),D));
end

n_dist = length(D);
n_targ = length(T);


k1 = numedges(G);
V_in_initial = []; % control on targets if those are directly connected to a disturbance
population = setdiff(setdiff(1:n, T), D);

for i = 1:n_targ % decouple targets from disturbances and other targets
    a = 0;
    [sout, ~] = findedge(G,inedges(G,T(i))); % find nodes with outgoing edges onto targets§
    if ~isempty(intersect(sout,D)) % check if those are disturbance nodes
        DT = intersect(sout,D);
        L1 = length(DT);
        V_in_initial = [V_in_initial, T(i)]; % if yes put the control on the target nodes
        G = rmedge(G,DT,T(i)); % remove incoming edges of targets
        a = 1;
    end
    if a
        [sout_2, ~] = findedge(G,inedges(G,T(i)));
        if isempty(sout_2)
            s = randsample(population, min(length(population),ceil(L1)));
            G = addedge(G, s, T(i), 1);
        end
    end
end

% k2 = numedges(G);
% 
% if abs(k2-k1) > 0
% 
%     s1 = randsample(population, min(length(population),abs(k2-k1))); % add as many edges as the removed ones
%     s2 = randsample(population, length(s1));
%     
%     while isequal(s1,s2)
%         s2 = randsample(population, min(length(population),length(s1)));
%     end
%     
%     for i = 1:length(s1)
%         G = addedge(G, s1(i), s2(i), 1);
%     end
% 
% end

% 
% k1 = numedges(G);
% V_in_initial = []; % control on targets if those are directly connected to a disturbance
% population = setdiff(setdiff(1:n, T), D);
% 
% for i = 1:n_targ % decouple targets from disturbances and other targets
%     a = 0;
%     [sout, ~] = findedge(G,inedges(G,T(i))); % find nodes with outgoing edges onto targets
%     if ~isempty(intersect(sout,D)) % check if those are disturbance nodes
%         DT = intersect(sout,D);
%         L1 = length(DT);
%         V_in_initial = [V_in_initial, T(i)]; % if yes put the control on the target nodes
%         %             G = rmedge(G,inedges(G,V_in_initial(end))); % remove incoming edges of targets where the input is directly applied
%         G = rmedge(G,DT,T(i)); % remove incoming edges of targets
%         a = 1;
%     end
%     if ~isempty(intersect(sout,T)) % check if those are target nodes
%         TT = intersect(sout,T);
%         L2 = length(TT);
%         G = rmedge(G,TT,T(i)); % remove incoming edges of targets
%         a = 1;
%     end
%     if a
%         [sout_2, ~] = findedge(G,inedges(G,T(i)));
%         if isempty(sout_2)
%             s = randsample(population, min(length(population),ceil((L1+L2)/2)));
%             G = addedge(G, s, T(i), 1);
%         end
%     end
% end
% 
% for i = 1:n_dist % decouple disturbances from other disturbances
%     a = 0;
%     [~,s_in] = findedge(G,outedges(G,D(i))); % find nodes with ingoing edges from disturbance
%     if ~isempty(intersect(s_in,D)) % check if those are disturbance nodes
%         DD = intersect(s_in,D);
%         L3 = length(DD);
%         G = rmedge(G,D(i),DD); % remove outcoming edges of disturbance
%         a = 1;
%     end
%     if a
%         [~,s_in_2] = findedge(G,outedges(G,D(i))); % find nodes with ingoing edges from disturbance
%         if isempty(s_in_2)
%             s = randsample(population, min(length(population),ceil((L1+L2)/2)));
%             G = addedge(G, D(i), s, 1);
%         end
%     end
% end
% 
% k2 = numedges(G);
% s1 = randsample(population, min(length(population),abs(k2-k1))); % add as many edges as the removed ones
% s2 = randsample(population, length(s1));
% 
% while isequal(s1,s2)
%     s2 = randsample(population, min(length(population),length(s1)));
% end
% 
% for i = 1:length(s1)
%     G = addedge(G, s1(i), s2(i), 1);
% end


figure; % After edge removal by action of V_in_initial on targets directly connected to disturbances

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

a = 1;

% clear; close all; clc;
% 
% A = zeros(17,17);
% D = [1 2];
% T = [13 14];
% A(3,1) = 1; A(3,2) = 1;
% A(4,3) = 1; A(5,3) = 1;
% A(6,4) = 1; A(6,5) = 1;
% A(7,6) = 1; A(8,6) = 1;
% A(9,6) = 1; A(10,7) = 1; A(11,8) = 1; A(11,9) = 1;
% A(12,10) = 1; A(12,11) = 1;
% A(13,12) = 1; A(14,12) = 1;
% A(6,15) = 1; A(15,3) = 1;
% A(16,6) = 1; A(12,16) = 1;
% A(17,6) = 1; A(12,17) = 1;
% G = digraph(A');
% figure;
% plot(G)
% 
% Gq = digraph(A);
% figure;
% plot(Gq)

% 
% % A = zeros(9,9);
% % A(2,1) = 1; A(3,1) = 1; A(4,2) = 1; A(4,3) = 1;
% % A(5,4) = 1; A(6,4) = 1; A(7,4) = 1; A(8,4) = 1;
% % A(9,5) = 1; A(9,6) = 1; A(9,7) = 1; A(9,8) = 1;
% % G = digraph(A');
% % D = 1;
% % T = 9;
% 
% figure;
% plot(G)


[V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair_backprop_new(G, D, T);



a = 1;

% [V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair(G, D, T);
% 
% a = 1;

k = 20; % select a (C,A,B)-pair
COM = compensator_ddp(G, D, T, V_out{k}, V_in{k}, S_Min{k}, Z_Max{k});

a = 1;

% % COM = compensator_ddp_old(G, D, T);
% 
% a = 1;
% 
% GG = minreal(tf(ss(COM.A_e,COM.D_e,COM.T_e,zeros(length(T),length(D)))));

syms t real
t = 1;
expo = round(COM.T_e*expm(sparse(COM.A_e)*t)*COM.D_e)

a = 1;

figure
plot(digraph(COM.A_e'))

a = 1;
% 
% [V_out, V_in, TotalCost] = cab_pair_old(G, D, T, 'all')
V_out = submincutDDSF_final2(G,D,T,'V_out');
[V_out_cell, G_cell] = mincutDDSF_all(G,D,T,V_out,'V_out','all');
% [V_out_cell6, G_cell6] = mincutDDSF_all(G,D,T,V_out,'V_in','all');
% [V_out_cell7, G_cell7] = mincutDDSF_all(G,D,T,V_out_cell,'V_in','all',G_cell{end});
[V_out_cell2, G_cell2] = mincutDDSF_all(G,D,T,V_out_cell,'V_out','all',G_cell{end});
[V_out_cell3, G_cell3] = mincutDDSF_all(G,D,T,[V_out_cell, V_out_cell2],'V_out','all',G_cell2{end});
[V_out_cell4, G_cell4] = mincutDDSF_all(G,D,T,[V_out_cell, V_out_cell2, V_out_cell3],'V_out','all',G_cell3{end});
[V_out_cell5, G_cell5] = mincutDDSF_all(G,D,T,[V_out_cell, V_out_cell2, V_out_cell3, V_out_cell4],'V_out','all',G_cell4{end});


COM = compensator_ddp(G, D, T);

a = 1;

GG = minreal(tf(ss(COM.A_e,COM.D_e,COM.T_e,zeros(length(T),length(D)))));

syms t real
t = 1;
expo = COM.T_e*expm(COM.A_e*t)*COM.D_e

a = 1;

figure
plot(digraph(COM.A_e'))

V_in = submincutDDSF_final2(G,D,T,'V_in');
V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');
[Vin_opt, Vout, C1] = constrained_optimal_solution(G,D,T,V_in_all,'V_out');
V_out = submincutDDSF_final2(G,D,T,'V_out');
V_out_all = mincutDDSF_all(G,D,T,V_out,'V_out','all');
[Vin, Vout_opt, C2] = constrained_optimal_solution(G,D,T,V_out_all,'V_in');
[V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(Vin_opt, Vout, C1, Vin, Vout_opt, C2);

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
