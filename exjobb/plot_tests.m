function plot_tests(tests, x_axis, options)
%PLOT_TEST Summary of this function goes here
%   Detailed explanation goes here
    arguments
        tests {mustBeCell}
        x_axis {mustBeVector}
        options.legend_entries {mustBeText} = "Null"
        options.legend_title {mustBeText} = "Null"
        options.graph_name {mustBeText} = "Null"
        options.x_label {mustBeText} = "Null"
        options.boxchart {mustBeNumericOrLogical} = false
    end

    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    for t = tests
        results_all(end+1,:) = t{1}.results_cost;
        results_time_all(end+1,:) = t{1}.results_time;
        results_trivial_all(end+1,:) = t{1}.results_trivial;
    end
    tests_plot(results_all, results_time_all, results_trivial_all, x_axis, options.legend_title, options.legend_entries, options.graph_name, options.x_label, options.boxchart)
end


function tests_plot(results_all, results_time_all, results_trivial_all, x_axis, legend_title, legend_entries, graph_name,  x_label, boxchart)

    mycolors = [
    238,64,53;
    243,155,54;
    123,192,67;
    3,146,207;
    17,0,255;
    175,56,255;
    ]./255;

    subplot(1,3,1);
    hold on;
    if boxchart
        for r = 1:size(results_all, 1)
            add2boxchart(results_all(r,:), string(x_axis(r)))
        end
    else
        mean_list = [];
        for r = 1:size(results_all, 1)
            mean_list(:,end+1) = mean(results_all(r,:)); 
        end
        plot(x_axis, mean_list, "-o", 'DisplayName', legend_entries)
        xlim([min(x_axis) max(x_axis)])
    end
    title("Cost")
    ylabel("Cost [-]")
    xlabel(x_label)
    ax = gca; 
    ax.ColorOrder = mycolors;
    hold off;

    subplot(1,3,2);
    hold on;
    if boxchart
        for r = 1:size(results_time_all, 1)
            add2boxchart(results_time_all(r,:), string(x_axis(r)))
        end
    else
        mean_list = [];
        for r = 1:size(results_time_all, 1)
            mean_list(:,end+1) = mean(results_time_all(r,:)); 
        end
        plot(x_axis, mean_list, "-o", 'DisplayName', legend_entries)
        xlim([min(x_axis) max(x_axis)])
    end
    title(graph_name + newline + "Runtime")
    ylabel("Time [s]")
    xlabel(x_label)
    ax = gca; 
    ax.ColorOrder = mycolors;
    hold off;
    
    subplot(1,3,3);
    hold on;
    if boxchart
        for r = 1:size(results_trivial_all, 1)
            add2boxchart(results_trivial_all(r,:), string(x_axis(r)))
        end
    else
        mean_list = [];
        for r = 1:size(results_trivial_all, 1)
            mean_list(:,end+1) = mean(results_trivial_all(r,:)); 
        end
        plot(x_axis, mean_list, "-o", 'DisplayName', legend_entries)
        xlim([min(x_axis) max(x_axis)])
        lgd = legend;
        title(lgd, legend_title)
    end
    title("Trivial solutions")
    ylabel("Index [-]")
    xlabel(x_label)
    ylim([0 1])
    ax = gca; 
    ax.ColorOrder = mycolors;
    hold off;
end

function add2boxchart(results, test_name)
    boxchart(categorical(1:numel(results), 1:numel(results), repmat(test_name, 1, numel(results))), results)
end