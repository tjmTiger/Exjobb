clear; clc; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load new_england_buss_data.mat;

% -----------
% Format data
% -----------
droop=0.05;
e=0.001;
J=1000;
f_0=60;
w_0=2*pi*f_0;
n=10;
m=length(K)-n;
E2=diag([e*ones(1,m) ones(1,n) J*w_0*ones(1,n)]);
c=1/droop*eye(n);
A2=[-K(1:m,1:m) -K(1:m,m+1:n+m) zeros(m,n); zeros(n,n+m) eye(n);  -K(m+1:n+m,1:end) -c ];
A2e=E2\A2;
G = digraph(A2e');

% G = G - diag(diag(G));  % remove self-loops (set diagonal to zeros)

% ---------------
% Defina Clusters
% ---------------
% 39 - 49 : generators (exclude 49 & 39!!!)
C1 = [48 47 40];
C2 = [41 42];
C3 = [46 45 43 44];

% C1 = [28 29 26 48 38 25 47 37 2 30 40];
% C2 = [5 7 8 6 41 31 11 12 13 10 32 42];
% C3 = [16 24 21 23 22 46 36 45 35 19 43 33 20 34 44];

% plot_clusters(G, C1, C2, C3)

% -----------
% Set D and T
% -----------
rng(1)
% T = randsample([C1 C2], 2)';
% D = randsample(C3, 2)';

results_cost_OF = zeros(1,60);
% results_trivial_OF = zeros(1,60);
results_cost_DF = zeros(1,60);
results_trivial_DF = zeros(1,60);
% All possible combinations with T in C1, C2 and D in C3
T_candidates = nchoosek([C1 C2],2);
D_candidates = nchoosek(C3,2);
i = 1;
for T = T_candidates'
    for D = D_candidates
        disp(i)
        % plot_system(G, T, D);
        [cost_OF, ~] = decouple_2(G, T, D, "output_feedback");
        [cost_DF, trivial_DF] = decouple_2(G, T, D, "dynamical_feedback", 'cost_OF', cost_OF);
        results_cost_OF(i) = cost_OF;
        results_cost_DF(i) = cost_DF;
        results_trivial_DF(i) = trivial_DF;
        i = i+1;
    end
end

figure;
subplot(1,2,1)
boxplot([results_cost_OF' results_cost_DF'], ["Output" "Dynamical"])
title("Cost")
ylabel("Cost")

subplot(1,2,2)
boxplot(results_trivial_DF, "Dynamical")
title("No Interaction Sol.")
ylabel("Index")


sgtitle("New England System")

fontsize(12,"points")
position = get(gcf, 'Position');
position = [100, 100, 600, 600];
% savefig(gcf, "figures_new/new_england")
% print(gcf, "figures_new/new_england" + ".eps", "-depsc")

% ---------
% Functions
% ---------
function [cost, triv] = decouple_2(G, T, D, ddp, options)
    arguments
        G
        T
        D
        ddp
        options.cost_OF
    end
    n_targ = numel(T);
    n_dist = numel(D);
    switch ddp
        case "state_feedback"
            t_start = tic;
            V_in = submincutDDSF_final2(G,D,T,'V_in');
            results_time = toc(t_start);
    
        case "output_feedback"
            t_start = tic;
            V_in = submincutDDSF_final2(G,D,T,'V_in');
            V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');
            [Vin_opt, Vout, C1] = constrained_optimal_solution(G,D,T,V_in_all,'V_out');
            V_out = submincutDDSF_final2(G,D,T,'V_out');
            V_out_all = mincutDDSF_all(G,D,T,V_out,'V_out','all');
            [Vin, Vout_opt, C2] = constrained_optimal_solution(G,D,T,V_out_all,'V_in');
            [V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(Vin_opt, Vout, C1, Vin, Vout_opt, C2);
            results_time = toc(t_start);

            if ~isempty(V_in_best) % if there is a solution, get one of those solutions
                V_in = V_in_best{1};
                V_out = V_out_best{1};
            else % if no solution, set V_in and V_out to empty
                V_in = [];
                V_out = [];
            end

            triv = 0;
            cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
            
    
        case "dynamical_feedback"
            numel_V_OF = options.cost_OF * ( n_targ + n_dist );
            t_start = tic;
            [V_out, V_in, ~, ~, ~] = cab_pair_backprop_new(G, D, T);
            results_time = toc(t_start);
    
            % How many of solutions with lower cost than OF have 0 trivial solutions?
            V_in_new = {};
            V_out_new = {};
            no_triv_count = 0;
            for i = 1:numel(V_in)
                if numel(V_in{i}) + numel(V_out{i}) <= numel_V_OF
                    V_in_new{end+1} = V_in{i};
                    V_out_new{end+1} = V_out{i};
                    if no_trivial_solutions(V_in{i}, V_out{i}, T, D)
                        no_triv_count = no_triv_count+1;
                        % disp("Found!!!!")
                    end
                end
            end
            
    
            if ~isempty(V_in) % if there is a solution, get smallest of those solutions
                V = cellfun(@numel, V_in) + cellfun(@numel, V_out);
                i = find(V==min(V));
    
                V_in = V_in{i};
                V_out = V_out{i};
            else % if no solution, set V_in and V_out to empty
                V_in = [];
                V_out = [];
            end
            
            if (numel(V_in_new) + numel(V_out_new)) ~= 0
                triv = 1 - (no_triv_count / (numel(V_in_new) + numel(V_out_new))); % calc_trivial_solutions([V_in V_out], [T; D]);
            else % Avoid division with 0
                triv = 0;
            end
            cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
    
        otherwise
            disp("ERROR: Invalid type of ddp. Possible options are: state_feedback, output_feedback or dynamical_feedback")
    end
end

function trivial_solutions = calc_trivial_solutions(V_in, T)
    [~,~,ic] = unique([V_in T']);
    a_counts = accumarray(ic,1);
    v_in_on_T = sum(a_counts(:,1)~=1);

    trivial_solutions = v_in_on_T/numel(V_in);
    if (v_in_on_T == 0) & (numel(V_in) == 0)
            trivial_solutions = 0;
    end
end

function out = no_trivial_solutions(V_in, V_out, T, D)
    out = calc_trivial_solutions([V_in V_out], [T; D]) == 0;
end

function plot_system(G, T, D)
    N = numnodes(G);

    figure; % After edge removal by action of V_in_initial on targets directly connected to disturbances
    p = plot(G,'b');
    title('$\mathcal{G}$')
    nodeColors = 1 * ones(N, 1);
    nodeColors(T) = 2;
    nodeColors(D) = 3;
    p.NodeCData = nodeColors;
    
    blue = [0 0 0.6]; green = [0 1 0]; red = [1 0 0];
    map = [blue; green; red];
    colormap(map)
    p.MarkerSize = 8;

    hold on;
    legendEntries = {'Other nodes', 'Disturbance', 'Target'};
    hBlue = scatter(nan, nan, 1, blue, 'filled');
    hRed = scatter(nan, nan, 1, red, 'filled');
    hGreen = scatter(nan, nan, 1, green, 'filled');
    legend([hBlue, hRed, hGreen], legendEntries, 'Location', 'best');
    hold off;
end


function plot_clusters(G, C1, C2, C3)
    N = numnodes(G);

    figure; % After edge removal by action of V_in_initial on targets directly connected to disturbances
    p = plot(G,'b');
    title('$\mathcal{G}$')
    nodeColors = 1 * ones(N, 1);
    nodeColors(C1) = 2;
    nodeColors(C2) = 3;
    nodeColors(C3) = 4;
    p.NodeCData = nodeColors;
    
    blue = [0 0 0.6]; green = [0 1 0]; red = [1 0 0]; yellow = [1 1 0];
    map = [blue; red; green; yellow];
    colormap(map)
    p.MarkerSize = 8;

    hold on;
    legendEntries = {'Other nodes', 'C1', 'C2', 'C3'};
    hBlue = scatter(nan, nan, 1, blue, 'filled');
    hRed = scatter(nan, nan, 1, red, 'filled');
    hGreen = scatter(nan, nan, 1, green, 'filled');
    hYellow = scatter(nan, nan, 1, yellow, 'filled');
    legend([hBlue, hRed, hGreen, hYellow], legendEntries, 'Location', 'best');
    hold off;
end