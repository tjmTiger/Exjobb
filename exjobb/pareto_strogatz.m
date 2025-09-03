function pareto_strogatz()
    disp("pareto_strogatz")
    % Ranges
    params1 = 1:3;                 % first argument
    params2 = linspace(0,1,101);  % second argument
    
    results = [];
    
    % loop through parameters
    for p1 = params1
        for p2 = params2
            [f1, f2] = DDP_watts_strogatz(p1, p2);  % returns [f1, f2]
            
            % [param1, param2, f1, f2]
            results = [results; p1, p2, f1, f2];
        end
    end
    % extract objectives
    obj1 = results(:,3);
    obj2 = results(:,4);
    
    % find pareto front finder
    is_pareto = true(size(obj1));
    for i = 1:length(obj1)
        for j = 1:length(obj1)
            if all([obj1(j), obj2(j)] <= [obj1(i), obj2(i)]) && ...
               any([obj1(j), obj2(j)] <  [obj1(i), obj2(i)])
                is_pareto(i) = false;
                break;
            end
        end
    end
    
    %% plot
    figure; hold on;
    
    % all samples in light blue
    scatter(obj1, obj2, 20, 'b', 'filled', 'MarkerFaceAlpha',0.3);
    
    % pareto points big dark blue
    paretoPoints = results(is_pareto, :);  % [mode, param, f1, f2]
    
    for k = 1:size(paretoPoints,1)
        p1 = paretoPoints(k,1);
        p2 = paretoPoints(k,2);
        f1 = paretoPoints(k,3);
        f2 = paretoPoints(k,4);
        color = [0 0.4470 0.7410];
        scatter(f1, f2, 60, color, 'filled');
        text(f1+0.01, f2, sprintf('%.2f; %.2f', p1, p2), 'FontSize',8, 'Color', color);
    end
    
    xlabel('Cost');
    ylabel('Trivial Sol.');
    legend('All Samples','Pareto','Location','best');
    title('Pareto Front for Watts Strogats');
    grid on;
    savefig(gcf, "figures_new/pareto_Strogatz")
    print(gcf, "figures_new/pareto_Strogatz" + ".eps", "-depsc")
end

function [cost, trivial] = DDP_watts_strogatz(k, beta)
    disp("k = " + k)
    disp("beta = " + beta)
    fract_targ = 0.1;
    fract_dist = 0.1;
    n = 100;
    params = num2cell([n, k, beta]);
    [costs, ~, trivials] = run_test( ...
        @watts_strogatz, ...
        params, ...
        "sample_size", 500, ...
        "fraction_targets", fract_targ, ...
        "fraction_disturbances", fract_dist, ...
        "ddp", "state_feedback", ...
        "seed", 1000 ...
        );
    cost = mean(costs);
    trivial = mean(trivials);
end

