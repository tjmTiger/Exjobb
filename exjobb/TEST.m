function TEST(test_function, variable1_name, variable1, variable2_name, variable2, options)
    arguments
        % Function defining the test.
        test_function

        % Will be plotted as seperate functions on same plot.
        variable1_name string
        variable1 {mustBeRow}

        % Will be plotted on x-axis
        variable2_name string
        variable2 {mustBeRow}
        
        % Tuning parameters
        options.sample_size {mustBeNumeric} = 100
        options.size {mustBeNumeric} = 100
        options.fract_targ double = 0.1
        options.fract_dist double = 0.1
        
        options.graph_generating_algorithm = ["Erdos Renyi"]% , "Watts Strogratz", "Scale Free"]
        options.ddp string  = "state_feedback"
        options.boxchart {mustBeNumericOrLogical} = false
        options.seed {mustBeNumeric} = 0
        options.default_p = [0.03, 2, 0.9, 1, 0, 0]; % p, k, pw, alpha, beta, gamma
    end


    for graph_generating_algorithm = {...
            {"Erdos Renyi", [options.size, options.default_p(1)]}, ...
            {"Watts Strogratz", [options.size, options.default_p(2), options.default_p(3)]}, ...
            {"Scale Free", [options.size, options.default_p(4), options.default_p(5), options.default_p(6), 1, 1]}}
        if ismember(graph_generating_algorithm{1}{1}, options.graph_generating_algorithm)
            results = {};
            old_results = {}; % for comparing DF with OF
            figure();
            for v1 = variable1
                old_results_count = 1; % for comparing DF with OF
                for v2 = variable2
                    if numel(old_results) >= old_results_count
                        options.old_results = old_results{old_results_count};
                    else
                        options.old_results = {};
                    end
                    results{end+1} = test_function(graph_generating_algorithm, v1, v2, options);
                end
                plot_results(results, variable2, "legend_title", erase(variable1_name, "_"), "legend_entries", erase(string(v1), "_"), "graph_name", graph_generating_algorithm{1}{1}, "x_label", erase(variable2_name, "_"), "boxchart", options.boxchart);
                old_results = results;
                results = {};
            end
            fontsize(12,"points")
            position = get(gcf, 'Position');
            position = [100, 100, 600, 600];
            savefig(gcf, "figures_new/" + variable1_name + "_" + variable2_name + "_" + erase(graph_generating_algorithm{1}{1}," "))
            print(gcf, "figures_new/" + variable1_name + "_" + variable2_name + "_" + erase(graph_generating_algorithm{1}{1}," ") + ".eps", "-depsc")
        end
    end
end






function plot_results(tests, x_axis, options)
%PLOT_TEST Summary of this function goes here
%   Detailed explanation goes here
    arguments
        tests cell
        x_axis {mustBeRow}
        options.legend_entries string
        options.legend_title string
        options.graph_name string
        options.x_label string
        options.boxchart {mustBeNumericOrLogical}
    end

    mycolors = [
    238,64,53;
    243,155,54;
    123,192,67;
    3,146,207;
    17,0,255;
    175,56,255;
    ]./255;

    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    for t = tests
        results_all(end+1,:) = t{1}.results_cost;
        results_time_all(end+1,:) = t{1}.results_time;
        results_trivial_all(end+1,:) = t{1}.results_trivial;
    end

    subplot(1,3,1);
    hold on;
    if options.boxchart
        for r = 1:size(results_all, 1)
            add2boxchart(results_all(r,:), string(x_axis(r)))
        end
    else
        mean_list = [];
        for r = 1:size(results_all, 1)
            mean_list(:,end+1) = mean(results_all(r,:)); 
        end
        plot(x_axis, mean_list, "-o", 'DisplayName', options.legend_entries)
        xlim([min(x_axis) max(x_axis)])
    end
    title("Cost")
    ylabel("Cost [-]")
    xlabel(options.x_label)
    ax = gca; 
    ax.ColorOrder = mycolors;
    hold off;

    subplot(1,3,2);
    hold on;
    if options.boxchart
        for r = 1:size(results_time_all, 1)
            add2boxchart(results_time_all(r,:), string(x_axis(r)))
        end
    else
        mean_list = [];
        for r = 1:size(results_time_all, 1)
            mean_list(:,end+1) = mean(results_time_all(r,:)); 
        end
        plot(x_axis, mean_list, "-o", 'DisplayName', options.legend_entries)
        xlim([min(x_axis) max(x_axis)])
    end
    title(options.graph_name + newline + "Runtime")
    ylabel("Time [s]")
    xlabel(options.x_label)
    ax = gca; 
    ax.ColorOrder = mycolors;
    hold off;
    
    subplot(1,3,3);
    hold on;
    if options.boxchart
        for r = 1:size(results_trivial_all, 1)
            add2boxchart(results_trivial_all(r,:), string(x_axis(r)))
        end
    else
        mean_list = [];
        for r = 1:size(results_trivial_all, 1)
            mean_list(:,end+1) = mean(results_trivial_all(r,:)); 
        end
        plot(x_axis, mean_list, "-o", 'DisplayName', options.legend_entries)
        xlim([min(x_axis) max(x_axis)])
        lgd = legend;
        title(lgd, options.legend_title)
    end
    title("Trivial solutions")
    ylabel("Index [-]")
    xlabel(options.x_label)
    ylim([0 1])
    ax = gca; 
    ax.ColorOrder = mycolors;
    hold off;
end


function add2boxchart(results, test_name)
    boxchart(categorical(1:numel(results), 1:numel(results), repmat(test_name, 1, numel(results))), results)
end