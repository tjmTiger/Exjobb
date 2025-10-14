function [results_cost, results_time, results_trivial] = run_test(algorithm, parameters, options)
% RUN_TEST
%   Description:
%       Run DDP on a sample size of random graphs created according to
%       algorithm provided with parameters passed in using
%       seed ∈ [seed, seed + sample_size].
%   Output Arguments:
%       results_cost    : Cost according to cost function defined in
%                         decouple() function.
%       results_time    : Time it took to decouple. Note: doesnt work
%                         properly when using multi threding
%       results_trivial : Ammount of trivial solutions.
%   Input Arguments:
%       algorithm       : @function for graph generating alorithm that
%                         returns a digraph object. (@sfg, @erods_renyi etc.)
%       parameters      : parameters to be passed into the alorithm func.
%       options         : seed - specify seed for rng
%                         sample_size - specify how many graphs each combination
%                                       of variables will be tested on.
%                         fraction_targets - what fraction of nodes will be chosen as targets
%                         fraction_disturbances - what fraction of nodes will be chosen as disturbances
%                         ddp - specify ddp algorithm default is state.
%                               Choices:
%                                 "state_feedback" (default)
%                                 "output_feedback"
%                                 "dynamical_feedback"
    arguments
        algorithm
        parameters {mustBeCell}
        options.sample_size {mustBeNumeric}
        options.seed {mustBeNumeric} = 1000 % NOTE: not the only place where default seed is chosen!!!!, change in TEST() too!!!
        options.fraction_targets single = 0.1
        options.fraction_disturbances single = 0.1
        options.ddp string
        options.old_results = {} % for comparing DF with OF
    end
    seed = options.seed;
    
    results_cost = zeros(1,options.sample_size);
    results_time = zeros(1,options.sample_size);
    results_trivial = zeros(1,options.sample_size);
    start_t = tic();
    parfor i = 1:options.sample_size % parallel processing for speed, change "parfor" to "for" when debugging (ty cant halt functions inside a parfor loop)
        % create randomly generated graph
        G = algorithm(parameters{:}, seed + i);
        % add targets and disturbances, decouple useing given ddp.
        [results_cost(i), results_time(i), results_trivial(i)] = decouple(G, options.fraction_targets, options.fraction_disturbances, "ddp", options.ddp, "seed", seed+i, "old_results", options.old_results, "parfor_i", i);
    end
    disp("run_test finished in: " + toc(start_t))
end