function result = test_different_fractions(graph_generating_algorithm, fract_dist, fract_targ, options)
% TEST_DIFFERENT_FRACTIONS
%   Description:
%       Tests how graphs are effected by fraction of targets and fraction
%       of disturbances. Fraction of distrubances in legend, fraction of targets on x-axis
%   Output Arguments:
%       result    : set containing cost, runtime and trivial solutions.
%   Input Arguments:
%       graph_generating_algorithm : "Erdos Renyi", "Watts Strogratz" or "Scale Free"
%       fract_dist: float, fraction of disturbances |D|/n
%       fract_targ: float, fraction of targets |T|/n
%       options   : sample_size, ddp, seed
arguments
    graph_generating_algorithm 
    fract_dist single
    fract_targ single
    options 
end
switch graph_generating_algorithm{1}{1}
    case "Erdos Renyi"
        graph_algorithm = @erdos_renyi;
    case "Watts Strogratz"
        graph_algorithm = @watts_strogatz;
    case "Scale Free"
        graph_algorithm = @sfg;
end
params = num2cell(graph_generating_algorithm{1}{2});
[result.results_cost, result.results_time, result.results_trivial] = run_test( ...
    graph_algorithm, ...
    params, ...
    "sample_size", options.sample_size, ...
    "fraction_targets", fract_targ, ...
    "fraction_disturbances", fract_dist, ...
    "ddp", options.ddp, ...
    "seed", options.seed ...
    );