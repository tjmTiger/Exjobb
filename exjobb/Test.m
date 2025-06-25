classdef Test < handle
%TEST Stores properties and results of a test with given parameters including cost,
%runtime and ammount of trivial solutions.
%   Detailed explanation goes here

properties
    % graphs properties
    graph_algorithm 
    graph_parameters {mustBeCell}
    graph_seed {mustBeNumeric} = 0

    % test properties
    number_of_tests {mustBeNumeric} = 200
    fraction_targets {mustBeFloat} = 0.1
    fraction_disturbances {mustBeFloat} = 0.1

    % other
    results_cost
    results_time
    results_trivial
end
methods
    function self = Test(graph_algorithm, graph_parameters, varargin)
        %TEST Construct an instance of this class
        %   Detailed explanation goes here

        self.graph_algorithm = graph_algorithm;
        self.graph_parameters = num2cell(graph_parameters);
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'number_of_tests'
                    self.number_of_tests = varargin{2};
                case 'fraction_targets'
                    self.fraction_targets = varargin{2};
                case 'fraction_disturbances'
                    self.fraction_disturbances = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
    end

    function run(self)
        %METHOD Perform a test according to object properties
        %   Detailed explanation goes here
        
        for i = 1:self.number_of_tests
            t_start = tic;
            G = self.graph_algorithm(self.graph_parameters{:},self.graph_seed + i);
            disp("new graph time: " + toc(t_start))
            t_start = tic;
            [self.results_cost(end+1), self.results_time(end+1), self.results_trivial(end+1)] = decouple(G, self.fraction_targets, self.fraction_disturbances);
            disp("decoupling time: " + toc(t_start))
        end
    end

    function add_to_plot(self)
        % pass
    end
end
end

%-----------------------------------------------%
%                                               %
%                 Local Functions               %
%                                               %
%-----------------------------------------------%

function [cost, runtime, trivial_solutions] = decouple(G, fract_targ, fract_dist, list_targ, list_dist)
    %DECOUPLE Add targets and ditrubances to a graph and run alorithm to
    %decouple them.
    %   Number of targets and disturbances is chosen according to frac_targ
    %   and frac_dist.
    %   Note: frac_targ + frac_dist <= 1
    
    if fract_targ + fract_dist > 1
        error('input argument invalid, fract_targ + fract_dist must be lesst than 1')
    end
    
    % cleaned up
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    
    N = numnodes(G);
    T = [];
    
    switch nargin
        case 3
            while isempty(T) % add targets and disturbances
                n_dist = ceil(fract_dist*N);
                n_targ = ceil(fract_targ*N);
                D = sort(randsample(N, n_dist));
                T = sort(randsample(setdiff(1:N', D), n_targ))';
            end

        case 5
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

        otherwise
            disp('input argument invalid')
    end
    
    n_dist = length(D);
    n_targ = length(T);
    
    V_in_initial = []; % control on targets if those are directly connected to a disturbance
    
    for i = 1:n_targ % decouple targets from disturbances and other targets
        [sout, ~] = findedge(G,inedges(G,T(i))); % find nodes with outgoing edges onto targets§
        if ~isempty(intersect(sout,D)) % check if those are disturbance nodes
            DT = intersect(sout,D);
            V_in_initial = [V_in_initial, T(i)]; % if yes put the control on the target nodes
            G = rmedge(G,DT,T(i)); % remove incoming edges of targets
        end
    end
    
    A = full(adjacency(G))';
    G = digraph(A');

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  State Feedback  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    t_start = tic;
    V_in = submincutDDSF_final2(G,D,T,'V_in');
    runtime = toc(t_start);
    
    [~,~,ic] = unique([V_in T']);
    a_counts = accumarray(ic,1);
    v_in_on_T = sum(a_counts(:,1)~=1);
    
    trivial_solutions = v_in_on_T/numel(V_in);
    if (v_in_on_T == 0) & (numel(V_in) == 0)
            trivial_solutions = 0;
    end
    
    cost = (2*numel(V_in)) / ( n_targ + n_dist );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Output Feedback  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % t_start = tic; % dubbelkolla vila av de nedan ska man mäta tid på.
    % 
    % V_in = submincutDDSF_final2(G,D,T,'V_in');
    % V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','all');
    % [Vin_opt, Vout, C1] = constrained_optimal_solution(G,D,T,V_in_all,'V_out');
    % V_out = submincutDDSF_final2(G,D,T,'V_out');
    % V_out_all = mincutDDSF_all(G,D,T,V_out,'V_out','all');
    % [Vin, Vout_opt, C2] = constrained_optimal_solution(G,D,T,V_out_all,'V_in');
    % [V_in_best, V_out_best, ~, ~] = global_constrained_optimal_solution(Vin_opt, Vout, C1, Vin, Vout_opt, C2);
    % 
    % runtime = toc(t_start);
    % V_in = cell2mat(V_in_best);
    % V_out = cell2mat(V_out_best);
    % 
    % [~,~,ic] = unique([V_in T']);
    % a_counts = accumarray(ic,1);
    % v_in_on_T = sum(a_counts(:,1)~=1);
    % 
    % trivial_solutions = v_in_on_T/numel(V_in);
    % if (v_in_on_T == 0) & (numel(V_in) == 0)
    %         trivial_solutions = 0;
    % end
    % 
    % cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Dynamical Feedback  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % cost = ( numel(V_in) + numel(V_out)) / ( n_targ + n_dist );
end