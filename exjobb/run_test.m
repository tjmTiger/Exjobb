function [results_cost, results_time, results_trivial] = run_test(algorithm, parameters, options)
%   Description:
%       Run DDP on a random graph created according to algorithm provided,
%       with parameters passed in.
%   Output Arguments:
%       results_cost    : Cost according to cost function defined in
%                         decouple() function.
%       results_time    : Time it took to decouple. Note: doesnt work
%                         properly when using multi threding
%       results_trivial : Ammount of trivial solutions.
%   Input Arguments:
%       algorithm       : @Function for graph generating alorithm that
%                         returns a digraph object.
%       parameters      : Parameters to be passed into the alorithm func.
%       options         : seed, sample_size, fraction_targets, fraction_disturbances

    arguments
        algorithm
        parameters {mustBeCell}
        options.sample_size {mustBeNumeric} = 200
        options.seed {mustBeNumeric} = 0
        options.fraction_targets {mustBeNumeric} = 0.1
        options.fraction_disturbances {mustBeNumeric} = 0.1
    end
    
    results_cost = zeros(1,options.sample_size);
    results_time = zeros(1,options.sample_size);
    results_trivial = zeros(1,options.sample_size);
    for i = 1:options.sample_size % parfor
        G = algorithm(parameters{:}, options.seed + i);
        [results_cost(i), results_time(i), results_trivial(i)] = decouple(G, options.fraction_targets, options.fraction_disturbances);
    end
end