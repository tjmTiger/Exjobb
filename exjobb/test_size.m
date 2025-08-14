function result = test_size(graph_generating_algorithm, ~, size, options)
    switch graph_generating_algorithm{1}{1}
        case "Erdos Renyi"
            graph_algorithm = @erdos_renyi;
        case "Watts Strogratz"
            graph_algorithm = @watts_strogatz;
        case "Scale Free"
            graph_algorithm = @sfg;
    end
    params = num2cell([size, graph_generating_algorithm{1}{2}(2:end)]);
    [result.results_cost, result.results_time, result.results_trivial] = run_test( ...
        graph_algorithm, ...
        params, ...
        "sample_size", options.sample_size, ...
        "fraction_targets", options.fract_targ, ...
        "fraction_disturbances", options.fract_dist, ...
        "ddp", options.ddp ...
        );



% 
% 
%     arguments
%         options.sample_size {mustBeNumeric} = 200
%         options.fract_targ {mustBeNumeric} = 0.3
%         options.fract_dist {mustBeNumeric} = 0.1
%         options.boxchart {mustBeNumericOrLogical} = true
%         options.ddp {mustBeText}  = "state_feedback"
%     end
%     for graph_generating_algorithm = {{"Erdos Renyi", [0.03]}, {"Watts Strogratz", [2, 0.2]}, {"Scale Free", [0.4, 0.2, 0.4, 1, 1]}}
%         % switch case needed for text used when plotting and saving
%         switch graph_generating_algorithm{1}{1}
%             case "Erdos Renyi"
%                 graph_algorithm = @erdos_renyi;
%             case "Watts Strogratz"
%                 graph_algorithm = @watts_strogatz;
%             case "Scale Free"
%                 graph_algorithm = @sfg;
%         end
%         tests = {};
%         figure();
%         for n = 30:30:180
%             params = num2cell([n, graph_generating_algorithm{1}{2}]);
%             [test.results_cost, test.results_time, test.results_trivial] = run_test( ...
%                 graph_algorithm, ...
%                 params, ...
%                 "sample_size", options.sample_size, ...
%                 "fraction_targets", options.fract_targ, ...
%                 "fraction_disturbances", options.fract_dist, ...
%                 "ddp", options.ddp ...
%                 );
%             tests{end+1} = test;
%         end
%         plot_tests(tests, 30:30:180, "graph_name", graph_generating_algorithm{1}{1}, "x_label", "Size", "boxchart", options.boxchart);
%         tests = {};
%         fontsize(12,"points")
%         position = get(gcf, 'Position');
%         position = [100, 100, 600, 600];
%         print(gcf, "figures_new/size_" + erase(graph_generating_algorithm{1}{1}," ") + ".eps", "-depsc")
%     end
% end
% 
