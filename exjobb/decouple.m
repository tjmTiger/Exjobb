function [cost, results_time, trivial_solutions] = decouple(G, fract_targ, fract_dist, options)
% DECOUPLE
%   Description:
%       Adds targets and disturbances, decouples using specified ddp.
%   Output Arguments:
%       results_cost    : Cost according to cost function depending on ddp
%       results_time    : Time it took to decouple.
%       results_trivial : Ammount of trivial solutions.
%   Input Arguments:
%       G               : digraph object
%       fract_targ      : float, fraction of targets
%                         numel(T) = nunodes(G)*fract_targ
%       fract_dist      : float, fraction of disturbances
%                         numel(D) = nunodes(G)*fract_dist
%                         Note: frac_targ + frac_dist <= 1
%       options         : ddp - specify ddp algorithm default is state.
%                               Choices:
%                                 "state_feedback" (default)
%                                 "output_feedback"
%                                 "dynamical_feedback"
%                         seed - specify seed for rng
%                         old_results, parfor_i - used for DF trivial solutions
arguments
    G
    fract_targ single
    fract_dist single
    options.ddp string = "state_feedback"
    options.seed {mustBeNumeric} = -1
    % options.list_targ = []
    % options.list_dist = []
    options.old_results
    options.parfor_i
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
        % no seed given
    otherwise
        rng(options.seed)
end

N = numnodes(G);
% choose targets and disturbances
T = [];
n_dist = ceil(fract_dist*N);
n_targ = ceil(fract_targ*N);
T = sort(randsample(N, n_targ));
D = sort(randsample(setdiff(1:N', T), n_dist))';
% plot_system(G, T, D)

n_dist = length(D);
n_targ = length(T);

A = full(adjacency(G))';
G = digraph(A');

% decouple using chosen ddp algorithm
switch options.ddp
    case "state_feedback"
        t_start = tic;
        V_in = submincutDDSF_final2(G,D,T,'V_in');

        % return:
        results_time = toc(t_start);
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
        
        if ~isempty(V_in_best) % if there is one or more solutions, pick one of those solutions
            V_in = V_in_best{1};
            V_out = V_out_best{1};
        else % if no solution, set V_in and V_out to empty
            V_in = [];
            V_out = [];
        end

        trivial_solutions = calc_trivial_solutions([V_in V_out], [T; D]);
        cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );

    case "dynamical_feedback"
        numel_V_OF = options.old_results.results_cost(options.parfor_i) * ( n_targ + n_dist ); % convert old cost back to number of V_in + V_out
        t_start = tic;
        [V_out, V_in, ~, ~, ~] = cab_pair_backprop_new(G, D, T);
        results_time = toc(t_start);

        % How many of solutions with lower cost than OF have 0 trivial
        % solutions? (index for DF)
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
            trivial_solutions = 1 - (no_triv_count / (numel(V_in_new) + numel(V_out_new))); % calc_trivial_solutions([V_in V_out], [T; D]);
        else % Avoid division with 0 (issue when V_in and V_out are empty)
            trivial_solutions = 0;
        end
        cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
    otherwise
        disp("ERROR: Invalid type of ddp. Possible options are: state_feedback, output_feedback or dynamical_feedback")
end
end

function trivial_solutions = calc_trivial_solutions(V, TD)
% calculates trivial solution according to index formula.
% for SF use calc_trivial_solutions(V_in, T)
% for OF use calc_trivial_solutions([V_in V_out], [T D])
    [~,~,ic] = unique([V TD']);
    a_counts = accumarray(ic,1);
    v_on_TD = sum(a_counts(:,1)~=1);

    trivial_solutions = v_on_TD/numel(V);
    if (v_on_TD == 0) & (numel(V) == 0)
            trivial_solutions = 0;
    end
end

function out = no_trivial_solutions(V_in, V_out, T, D)
% return true if (V_in U V_out) ∩ (T U D) = ∅, else return false
    out = calc_trivial_solutions([V_in V_out], [T; D]) == 0;
end

function plot_system(G, T, D)
% Plot with targets and disturbances marked
% G: digraph object
% T: array of targets
% D: array of disturbances
    N = numnodes(G);

    figure;
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
    legendEntries = {'Disturbance', 'Target', 'Other nodes'};
    hRed = scatter(nan, nan, 100, red, 'filled');
    hGreen = scatter(nan, nan, 100, green, 'filled');
    hYellow = scatter(nan, nan, 100, blue, 'filled');
    legend([hRed, hGreen, hYellow], legendEntries, 'Location', 'best');
    hold off;
end