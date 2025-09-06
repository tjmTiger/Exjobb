function [cost, results_time, trivial_solutions] = decouple(G, fract_targ, fract_dist, options)
% Note: frac_targ + frac_dist <= 1
% Choices for ddp: 
%   state_feedback
%   output_feedback
%   dynamical_feedback
arguments
    G 
    fract_targ 
    fract_dist
    options.ddp {mustBeText} = "state_feedback"
    options.seed {mustBeNumeric} = -1
    options.list_targ = []
    options.list_dist = []
end

% cleaned up
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

if fract_targ + fract_dist > 1
    error('Invalid argument list. fract_targ + fract_dist must be lesst than 1')
end

switch options.seed
    case -1
        % pass
    otherwise
        rng(options.seed)
end

N = numnodes(G);
T = [];
if isempty(options.list_targ)
    while isempty(T) % add targets and disturbances
        n_dist = ceil(fract_dist*N);
        n_targ = ceil(fract_targ*N);
        T = sort(randsample(N, n_dist));
        D = sort(randsample(setdiff(1:N', T), n_targ))';
        T = sort(setdiff(randsample(N, n_targ),D));
    end
else
    if isempty(options.list_dist)
        disp("ERROR, Invalid input, list_dist cant be empty if list_targ is not empty.")
        return
    end
    list_dist = options.list_dist;
    list_targ = options.list_targ;
    % if length(list_targ) > length(list_dist)
    %     list_targ = setdiff(list_targ, list_dist);
    % else
    %     list_dist = setdiff(list_dist, list_targ);
    % end

    while isempty(T) % add targets and disturbances
        n_dist = ceil(fract_dist*N);
        n_targ = ceil(fract_targ*N);
        if n_dist > length(list_dist)
            disp("WARNING: Request " + n_dist + " distubances, but only " + length(list_dist) + " provided!")
            n_dist = length(list_dist);
        elseif n_targ > length(list_targ)
            disp("WARNING: Request " + n_targ + " targets, but only " + length(list_targ) + " provided!")
            n_targ = length(list_targ);
        end
        D = sort(randsample(list_dist, n_dist));
        T = sort(setdiff(randsample(N, n_targ),D));
    end
end

n_dist = length(D);
n_targ = length(T);

V_in_initial = []; % control on targets if those are directly connected to a disturbance
population = setdiff(setdiff(1:N, T), D);

% figure; % After edge removal by action of V_in_initial on targets directly connected to disturbances
% p = plot(G,'b');
% title('$\mathcal{G}$')
% nodeColors = 1 * ones(N, 1); % Default to value 3 (Yellow)
% nodeColors(T) = 2;
% nodeColors(D) = 3;
% p.NodeCData = nodeColors;
% colormap(jet); % Use the 'jet' colormap
% p.MarkerSize = 8; % Increase or decrease the size of the nodes
% hold on; % Hold on to the current plot
% legendEntries = {'Disturbance', 'Target', 'Other nodes'};
% hRed = scatter(nan, nan, 100, 'r', 'filled'); % Placeholder for red nodes
% hGreen = scatter(nan, nan, 100, 'g', 'filled'); % Placeholder for green nodes
% hYellow = scatter(nan, nan, 100, 'b', 'filled'); % Placeholder for yellow nodes
% legend([hRed, hGreen, hYellow], legendEntries, 'Location', 'best');
% hold off; % Release the hold on the current plot

A = full(adjacency(G))';
G = digraph(A');

switch options.ddp
    case "state_feedback"
        t_start = tic;
        V_in = submincutDDSF_final2(G,D,T,'V_in');
        results_time = toc(t_start);
        
        % V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');
        % check and display how many control nodes are placed on target nodes
        % v_in_on_T = numel(V_in);
        % for v_in = V_in_all
        %     v_in = v_in{1};
        %     [~,~,ic] = unique([v_in T']);
        %     a_counts = accumarray(ic,1);
        %     v_in_on_T_next = sum(a_counts(:,1)~=1);
        %     if v_in_on_T > v_in_on_T_next
        %             v_in_on_T = v_in_on_T_next;
        %     end
        % end
        
        trivial_solutions = calc_trivial_solutions(V_in, T);
        cost = (2*numel(V_in)) / ( n_targ + n_dist );

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

        V_in = cell2mat(V_in_best);
        V_out = cell2mat(V_out_best);

        trivial_solutions = calc_trivial_solutions(V_in, T);
        cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );

    case "dynamical_feedback"
        t_start = tic;
        [V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair_backprop_new(G, D, T);
        results_time = toc(t_start);
        V_in = cell2mat(V_in)';
        V_out = cell2mat(V_out)';
        trivial_solutions = calc_trivial_solutions(V_in, T);
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