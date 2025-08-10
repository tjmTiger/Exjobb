function test_strogatz_fraction_beta(options)
    arguments
        options.sample_size {mustBeNumeric} = 200
        options.n {mustBeNumeric} = 100
    end
    [options.n, 2, 0.2]
    graph_algorithm = @watts_strogatz;
    tests = {};
    figure();
    params = num2cell(graph_generating_algorithm{1}{2});
    for fraction = 0.05:0.05:0.3
        for beta = 0.05:0.05:0.3
            [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
                graph_algorithm, ...
                params, ...
                "sample_size", options.sample_size, ...
                "fraction_targets", fraction, ...
                "fraction_disturbances", fraction ...
                );
            tests{end+1} = test;
        end
        plot_tests(tests, 0.05:0.05:0.3,"legend_title", "Dist.Frac", "legend_entries", string(fract_targ), "graph_name", "Watts Strogratz", "x_label", "Targ.Frac.");
        tests = {};
    end
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    print(gcf, "figures_new/fraction_beta_" + erase("Watts Strogratz"," ") + ".eps", "-depsc")
end

