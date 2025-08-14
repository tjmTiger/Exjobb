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
    options.list_targ = 'Null'
    options.list_dist = 'Null'
end

if fract_targ + fract_dist > 1
    error('Invalid argument list. fract_targ + fract_dist must be lesst than 1')
end

% cleaned up
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

N = numnodes(G);

% fract_targ = 0.001;
% fract_dist = 0.001;
T = [];

switch options.list_targ
    case 'Null'
        while isempty(T) % add targets and disturbances
            n_dist = ceil(fract_dist*N);
            n_targ = ceil(fract_targ*N);
            D = sort(randsample(N, n_dist));
            T = sort(randsample(setdiff(1:N', D), n_targ))';
            % T = sort(setdiff(randsample(N, n_targ),D));
        end
    otherwise
        if options.list_dist == 'Null'
            disp("ERROR, Invalid input, list_dist cant be null if list_targ is not null.")
            return
        end
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
        
        [~,~,ic] = unique([V_in T']);
        a_counts = accumarray(ic,1);
        v_in_on_T = sum(a_counts(:,1)~=1);
        
        trivial_solutions = v_in_on_T/numel(V_in);
        if (v_in_on_T == 0) & (numel(V_in) == 0)
                trivial_solutions = 0;
        end
        
        cost = (2*numel(V_in)) / ( n_targ + n_dist );

    case "output_feedback"
        t_start = tic; % dubbelkolla vila av de nedan ska man mäta tid på.

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

        [~,~,ic] = unique([V_in T']);
        a_counts = accumarray(ic,1);
        v_in_on_T = sum(a_counts(:,1)~=1);

        trivial_solutions = v_in_on_T/numel(V_in);
        if (v_in_on_T == 0) & (numel(V_in) == 0)
                trivial_solutions = 0;
        end
        cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
    case "dynamical_feedback"
        % todo: write code for dynamical feedback
        % cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
    otherwise
        disp("ERROR: Invalid type of ddp. Possible options are: state_feedback, output_feedback or dynamical_feedback")
end