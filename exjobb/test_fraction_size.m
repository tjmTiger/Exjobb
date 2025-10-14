function result = test_fraction_size(graph_generating_algorithm, n, fraction, options)
% TEST_FRACTION_SIZE
%   Description:
%       Tests how differently sized graphs are effected by fraction of
%       targets/disturbances. Size (n) in legend, fraction on x-axis
%   Output Arguments:
%       result    : set containing cost, runtime and trivial solutions.
%   Input Arguments:
%       graph_generating_algorithm : "Erdos Renyi", "Watts Strogratz" or "Scale Free"
%       n         : integer, graph size, number of nodes, |V|
%       fractions : float, fraction of targets and disturbances |T|/n, |D|/n
%       options   : sample_size, ddp, seed
arguments
    graph_generating_algorithm 
    n {mustBeNumeric}
    fraction single
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
params = num2cell([n, graph_generating_algorithm{1}{2}(2:end)]);
[result.results_cost, result.results_time, result.results_trivial] = run_test( ...
    graph_algorithm, ...
    params, ...
    "sample_size", options.sample_size, ...
    "fraction_targets", fraction, ...
    "fraction_disturbances", fraction, ...
    "ddp", options.ddp, ...
    "seed", options.seed ...
    );