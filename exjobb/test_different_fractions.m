function test_different_fractions(options)
    arguments
        options.sample_size {mustBeNumeric} = 200
        options.n {mustBeNumeric} = 100
    end
    for graph_generating_algorithm = {{"Erdos Renyi", [options.n, 0.03]}, {"Watts Strogratz", [options.n, 2, 0.2]}, {"Scale Free", [options.n, 0.4, 0.2, 0.4, 1, 1]}}
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
        params = num2cell(graph_generating_algorithm{1}{2});
        for fract_dist = 0.05:0.05:0.3
            for fract_targ = 0.05:0.05:0.3
                [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
                    graph_algorithm, ...
                    params, ...
                    "sample_size", options.sample_size, ...
                    "fraction_targets", fract_targ, ...
                    "fraction_disturbances", fract_dist ...
                    );
                tests{end+1} = test;
            end
            plot_tests(tests, 0.05:0.05:0.3,"legend_title", "Dist.Frac", "legend_entries", string(fract_targ), "graph_name", graph_generating_algorithm{1}{1}, "x_label", "Targ.Frac.");
            tests = {};
        end
        fontsize(12,"points")
        position = get(gcf, 'Position');
        position = [100, 100, 600, 600];
        print(gcf, "figures_new/different_fractions_" + erase(graph_generating_algorithm{1}{1}," ") + ".eps", "-depsc")
    end
end

