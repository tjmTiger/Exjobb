function test_node_degree(options)
    arguments
        options.sample_size {mustBeNumeric} = 200
        options.n {mustBeNumeric} = 100
        options.fract_targ {mustBeNumeric} = 0.3
        options.fract_dist {mustBeNumeric} = 0.1
        options.boxchart {mustBeNumericOrLogical} = true
    end
    
    % Erdos Renyi
    graph_algorithm = @erdos_renyi;
    tests = {};
    figure();
    for i = 2:1:7
        params = num2cell([options.n, 2.2*i/(options.n-1)]);
        [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
            graph_algorithm, ...
            params, ...
            "sample_size", options.sample_size, ...
            "fraction_targets", options.fract_targ, ...
            "fraction_disturbances", options.fract_dist ...
            );
        tests{end+1} = test;
    end
    plot_tests(tests, 2:1:7, "graph_name", "Erdos Renyi", "x_label", "Average degree", "boxchart", options.boxchart);
    tests = {};
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    print(gcf, "figures_new/node_degree_" + erase("Erdos Renyi"," ") + ".eps", "-depsc")

    % Watts Strogratz
    graph_algorithm = @watts_strogatz;
    tests = {};
    figure();
    for i = 2:2:7
        params = num2cell([options.n, floor(i/2), 0.2]);
        [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
            graph_algorithm, ...
            params, ...
            "sample_size", options.sample_size, ...
            "fraction_targets", options.fract_targ, ...
            "fraction_disturbances", options.fract_dist ...
            );
        tests{end+1} = test;
    end
    plot_tests(tests, 2:2:7, "graph_name", "Watts Strogratz", "x_label", "Average degree", "boxchart", options.boxchart);
    tests = {};
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    print(gcf, "figures_new/node_degree_" + erase("Watts Strogratz"," ") + ".eps", "-depsc")

    % Scale Free
    graph_algorithm = @sfg;
    tests = {};
    figure();
    for beta = [0.55 0.70 0.79 0.845 0.871 0.893]
        alpha = (1-beta)/2;
        gamma = (1-beta)/2;
        delta_in = 10;
        delta_out = 10;
        params = num2cell([options.n, alpha, beta, gamma, delta_in, delta_out]);
        [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
            graph_algorithm, ...
            params, ...
            "sample_size", options.sample_size, ...
            "fraction_targets", options.fract_targ, ...
            "fraction_disturbances", options.fract_dist ...
            );
        tests{end+1} = test;
    end
    plot_tests(tests, 2:1:7, "graph_name", "Scale Free", "x_label", "Average degree", "boxchart", options.boxchart);
    tests = {};
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    print(gcf, "figures_new/node_degree_" + erase("Scale Free"," ") + ".eps", "-depsc")
end

