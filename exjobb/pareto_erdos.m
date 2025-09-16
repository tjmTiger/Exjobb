function pareto_erdos(ddp, n)
    disp("pareto_erdos")
    % Ranges
    % params1 = linspace(0.03,1,98);  % first argument
    params1 = logspace(0, -1.5229, 100);
    
    results = [];
    
    % loop through parameters
    for p1 = params1
        [f1, f2] = DDP_erods_renyi(p1, ddp, n);  % returns [f1, f2]
        
        % [param1, param2, f1, f2]
        results = [results; p1, 0, f1, f2];
    end
    % extract objectives
    obj1 = results(:,3);
    obj2 = results(:,4);
    
    % find pareto front
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
        fprintf('%.2f; %.2f\n', p1, p2)
    end
    
    xlabel('Cost');
    ylabel('Trivial Sol.');
    legend('All Samples','Pareto','Location','best');
    title('Pareto Front for Erdos Renyi');
    grid on;
    savefig(gcf, "figures_new/pareto_Erdos_" + ddp)
    print(gcf, "figures_new/pareto_Erdos_" + ddp + ".eps", "-depsc")
end

function [cost, triv] = DDP_erods_renyi(p, ddp, n)
    disp("p = " + p)
    fract_targ = 0.1;
    fract_dist = 0.1;
    params = num2cell([n, p]);
    [costs, ~, triv] = run_test( ...
        @erdos_renyi, ...
        params, ...
        "sample_size", 500, ...
        "fraction_targets", fract_targ, ...
        "fraction_disturbances", fract_dist, ...
        "ddp", ddp, ...
        "seed", 1000 ...
        );
    cost = mean(costs);
    triv = mean(triv);
end