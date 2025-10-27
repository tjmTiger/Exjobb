function result = test_size(graph_generating_algorithm, ~, n, options)
% TEST_SIZE
%   Description:
%       Tests how graphs are effected by their sizes. Note: Erdos Renyi p = 0.1 regardless of default params
%       due to graph size 180 > n > 30
%   Output Arguments:
%       result    : set containing cost, runtime and trivial solutions.
%   Input Arguments:
%       graph_generating_algorithm : "Erdos Renyi", "Watts Strogratz" or "Scale Free"
%       ~         : ignored input
%       n         : integer, integer, graph size, number of nodes, |V|
%       options   : sample_size, fract_targ, fract_dist, ddp, seed
arguments
    graph_generating_algorithm 
    ~
    n 
    options 
end
switch graph_generating_algorithm{1}{1}
    case "Erdos Renyi"
        graph_algorithm = @erdos_renyi;
        graph_generating_algorithm{1}{2}(2) = 0.1; % use p = 0.1 instead of 0.03
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
    "fraction_targets", options.fract_targ, ...
    "fraction_disturbances", options.fract_dist, ...
    "ddp", options.ddp, ...
    "seed", options.seed ...
    );