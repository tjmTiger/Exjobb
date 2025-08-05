clear; clc;

% Use Parallel pools to speed up the code.
p = gcp("nocreate");
NumWorkers = 8;
if isempty(p) % If no pool, create new one.
    parpool(NumWorkers)
else % If there is a pool, but its small, delete it and create new one.
    if p.NumWorkers < NumWorkers
        delete(gcp('nocreate'))
        parpool(NumWorkers)
    end
end

%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

% fraction and size
for graph_generating_algorithm = {{"Watts Strogratz", [2, 0.2]}, {"Scale Free", [0.4, 0.2, 0.4, 1, 1]}} % {{"Erdos Renyi", [0.03]}, {"Watts Strogratz", [2, 0.2]}, {"Scale Free", [0.4, 0.2, 0.4, 1, 1]}}
    switch graph_generating_algorithm{1}{1}
        case "Erdos Renyi"
            graph_algorithm = @erdos_renyi;
        case "Watts Strogratz"
            graph_algorithm = @watts_strogatz;
        case "Scale Free"
            graph_algorithm = @sfg;
    end
    tests = {};
    figure();
    for n = 100:50:200
        params = num2cell([n, graph_generating_algorithm{1}{2}]);
        for fraction = 0.05:0.05:0.3
            [results_cost, results_time, results_trivial] = run_test( ...
                graph_algorithm, ...
                params, ...
                "sample_size", 200, ...
                "fraction_targets", fraction, ...
                "fraction_disturbances", fraction ...
                );
            test.results_cost = results_cost;
            test.results_time = results_time;
            test.results_trivial = results_trivial;
            tests{end+1} = test;
        end
        plot_tests(tests, 0.05:0.05:0.3, "display_name", string(n), "graph_name", graph_generating_algorithm{1}{1});
        tests = {};
    end
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    print(gcf, "figures_new/size_fractions_" + erase(graph_generating_algorithm{1}{1}," ") + ".eps", "-depsc")
end
