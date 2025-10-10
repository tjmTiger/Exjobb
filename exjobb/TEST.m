function TEST(test_function, variable1_name, variable1, variable2_name, variable2, options)
% TEST
%   Description:
%       Runs given test_function in a tripple for loop:
%       for g = ["Erdos Renyi" , "Watts Strogratz", "Scale Free"]
%           for v1 = variable1
%               for v2 = variable2
%                   test_function(g, v1, v2, options)
%               end
%           end
%       end
%       NOTE: requires folder named "figures_new" in same directory to save
%             figures to. (otherwise comment out rows with savefig() and print()
%   Output Arguments:
%       None (plots results and saves figures)
%   Input Arguments:
%       test_function   : function defining what the test will do
%       variable1_name  : string, visual only for figure plot and .fig file name
%       variable1       : row array, will be plotted as separate lines, can be set to
%                         empty string if only one varable is to be plotted.
%       variable2_name  : string, visual only for figure plot and .fig file name
%       variable2       : row array, will be plotted on x-axis
%       options         : sample_size, size, fract_targ, fract_dist,
%                         graph_generating_algorithm - "Erdos Renyi", "Watts Strogratz" or "Scale Free"
%                         ddp - specify ddp algorithm default is state.
%                               Choices:
%                                 "state_feedback" (default)
%                                 "output_feedback"
%                                 "dynamical_feedback"
%                         boxchart - bool
%                                    true: present data in a box and wiskers diagram.
%                                    false: present data in normal plot. (default)
%                         seed - specify seed for rng
%                         default_p - array of default parameters [p, k, p_ws, alpha, beta, gamma]
    arguments
        % Function defining the test.
        test_function

        % Will be plotted as seperate functions on same plot.
        variable1_name string % only visual for figure plot and .fig file name
        variable1 {mustBeRow}

        % Will be plotted on x-axis
        variable2_name string % only visual for figure plot and .fig file name
        variable2 {mustBeRow}
        
        % Tuning parameters
        options.sample_size {mustBeNumeric} = 100
        options.size {mustBeNumeric} = 100
        options.fract_targ single = 0.1
        options.fract_dist single = 0.1
        
        options.graph_generating_algorithm {mustBeText} = ["Erdos Renyi" , "Watts Strogratz", "Scale Free"]
        options.ddp string  = "state_feedback"
        options.boxchart {mustBeNumericOrLogical} = false
        options.seed {mustBeNumeric} = 1000
        % default [p, k, pw, alpha, beta, gamma] parameters for random graph generation, note: not always used
        options.default_p {mustBeNumeric} = [0.03, 2, 0.9, 1, 0, 0]; % chosen according to pareto, [p, k, p_ws, alhpa. beta, gamma]
    end


    for graph_generating_algorithm = {...
            {"Erdos Renyi", [options.size, options.default_p(1)]}, ...
            {"Watts Strogratz", [options.size, options.default_p(2), options.default_p(3)]}, ...
            {"Scale Free", [options.size, options.default_p(4), options.default_p(5), options.default_p(6), 1, 1]}}
        if ismember(graph_generating_algorithm{1}{1}, options.graph_generating_algorithm)
            results = {};
            old_results = {}; % for comparing DF with OF, contains all results from previous variable1 loop
            figure();
            for v1 = variable1 % variable1 is plotted as separate lines
                old_results_count = 1; % for comparing DF with OF
                for v2 = variable2 % variable2 is plotted on x-axis
                    if numel(old_results) >= old_results_count % pass correct OF test result to DF test (from same graph)
                        options.old_results = old_results{old_results_count};
                        old_results_count = old_results_count+1;
                    else
                        options.old_results = {};
                    end
                    % perform test according to given test_function.
                    results{end+1} = test_function(graph_generating_algorithm, v1, v2, options);
                end
                plot_results(results, variable2, "legend_title", variable1_name, "legend_entries", string(v1), "graph_name", graph_generating_algorithm{1}{1}, "x_label", variable2_name, "boxchart", options.boxchart);
                old_results = results;
                results = {};
            end
            position = get(gcf, 'Position');
            position = [100, 100, 600, 600];
            % save figures, comment out or add folder named "figures_new" to direcotry with this function
            illegal_chars = [" ","\","}","{","^","$"];
            savefig(gcf, "figures_new/" + erase(options.ddp, illegal_chars) + "_" + erase(variable1_name, illegal_chars) + "_" + erase(variable2_name, illegal_chars) + "_" + erase(graph_generating_algorithm{1}{1},illegal_chars))
        end
    end
end



function plot_results(tests, x_axis, options)
% PLOT_RESULTS Visualization of TEST results.
%   test: y-values for the plot, set containing results with costs, runtimes and trivial solutions
%   x_axis: x-values for the plot, varable2
    arguments
        tests cell
        x_axis {mustBeRow}
        options.legend_entries string
        options.legend_title string
        options.graph_name string
        options.x_label string
        options.boxchart {mustBeNumericOrLogical}
    end
    
    % custom colors for plotted lines
    mycolors = [
    238,64,53;
    243,155,54;
    123,192,67;
    3,146,207;
    17,0,255;
    175,56,255;
    ]./255;
    
    % unpack results
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    for t = tests
        results_all(end+1,:) = t{1}.results_cost;
        results_time_all(end+1,:) = t{1}.results_time;
        results_trivial_all(end+1,:) = t{1}.results_trivial;
    end
    
    % Cost
    subplot(1,3,1);
    hold on;
    if options.boxchart % if boxchart option selected plot box and wisker diagram, else do normal plot
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
    
    % Runtime
    subplot(1,3,2);
    hold on;
    if options.boxchart % if boxchart option selected plot box and wisker diagram, else do normal plot
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
    
    % Trivial solutions
    subplot(1,3,3);
    hold on;
    if options.boxchart % if boxchart option selected plot box and wisker diagram, else do normal plot
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

    interpreter = 'latex';
    all_text = findall(gcf, "-property", "Interpreter");
    set(all_text, "Interpreter", interpreter)

    all_text = findall(gcf, "-property", "TickLabelInterpreter");
    set(all_text, "TickLabelInterpreter", interpreter)

    all_text = findall(gcf, "-property", "Fontsize");
    set(all_text, "Fontsize", 12)
end


function add2boxchart(results, test_name) % used in plot_results if boxchart = true
    boxchart(categorical(1:numel(results), 1:numel(results), repmat(test_name, 1, numel(results))), results)
end