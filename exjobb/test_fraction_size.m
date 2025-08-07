function test_fraction_size(options)
    arguments
        options.sample_size {mustBeNumeric} = 200
    end
    for graph_generating_algorithm = {{"Erdos Renyi", [0.03]}, {"Watts Strogratz", [2, 0.2]}, {"Scale Free", [0.4, 0.2, 0.4, 1, 1]}}
        % switch case needed for text used when plotting and saving
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
                [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
                    graph_algorithm, ...
                    params, ...
                    "sample_size", options.sample_size, ...
                    "fraction_targets", fraction, ...
                    "fraction_disturbances", fraction ...
                    );
                tests{end+1} = test;
            end
            plot_tests(tests, 0.05:0.05:0.3, "legend_title", "Size", "legend_entries", string(n), "graph_name", graph_generating_algorithm{1}{1}, "x_label", "Fractions");
            tests = {};
        end
        fontsize(12,"points")
        position = get(gcf, 'Position');
        position = [100, 100, 600, 600];
        print(gcf, "figures_new/size_fractions_" + erase(graph_generating_algorithm{1}{1}," ") + ".eps", "-depsc")
    end
end

