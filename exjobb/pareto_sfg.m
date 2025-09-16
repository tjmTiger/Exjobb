function pareto_sfg(ddp, n)
    disp("pareto_sfg")
    % Ranges
    % state
    params1 = linspace(0.1,1,10);  % first argument (alpha)
    params2 = linspace(0,1,11);  % second argument (beta)
    % params3 = linspace(0,1,11); % third argument (gamma)
    
    results = [];
    params_used = [];
    % loop through parameters
    for p1 = params1
        for p2 = params2
            if p1 + p2 <= 1
                params_used(end+1,:) = [p1,p2];
                [f1, f2] = DDP_sfg(p1, p2, ddp, n);  % returns [f1, f2]
                
                % [param1, param2, f1, f2]
                results = [results; p1, p2, f1, f2];
                % disp(f1 + ", " + f2)
            end
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
    figure;
    
    % all samples in light blue
    % scatter(obj1, obj2, 20, 'b', 'filled', 'MarkerFaceAlpha',0.3);
    % plot(obj1, obj2, "k-")
    colors = [255 0 0; 255 125 0; 255 255 0; 125 255 0; 0 255 0; 0 255 125; 0 255 255; 0 125 255; 0 0 255; 125 0 255; 255 0 255; 255 0 125]./255;

    subplot(3,1,1)
    hold on
    title("alpha")
    grid on;
    p1_last = params_used(1,1);
    obj = [];
    c = 1;
    for i = 1:numel(params_used(:,1))
        p1 = params_used(i,1);
        if p1 ~= p1_last
            color = colors(c,:);
            plot(obj(:,1), obj(:,2), "o-",'Color',color,'MarkerSize',4, 'MarkerFaceColor', color)
            obj = [];
            c = min(c+1, 12);
            p1_last = p1;
        end
        obj = [obj; obj1(i) obj2(i)];
    end
    % pareto points big dark blue
    paretoPoints = results(is_pareto, :);  % [mode, param, f1, f2]

    for k = 1:size(paretoPoints,1)
        p1 = paretoPoints(k,1);
        p2 = paretoPoints(k,2);
        f1 = paretoPoints(k,3);
        f2 = paretoPoints(k,4);
        color = [0 0.4470 0.7410];
        scatter(f1, f2, 60, color, 'filled');
        text(f1+0.001, f2, sprintf('%.2f; %.2f; %.2f\n', p1, p2, 1-p1-p2), 'FontSize',8, 'Color', color);
        fprintf('%.2f; %.2f; %.2f\n', p1, p2, 1-p1-p2)
    end
    hold off

    subplot(3,1,2)
    hold on
    title("beta")
    grid on;
    obj = [];
    c = 1;
    for p2_last = params2
        for i = 1:numel(params_used(:,1))
            p2 = params_used(i,2);
            if p2 == p2_last
                obj = [obj; obj1(i) obj2(i)];
            end
        end
        color = colors(c,:);
        if ~isempty(obj)
            plot(obj(:,1), obj(:,2), "o-",'Color',color,'MarkerSize',4, 'MarkerFaceColor', color)
        end
        obj = [];
        c = min(c+1, 12);
    end
    % pareto points big dark blue
    paretoPoints = results(is_pareto, :);  % [mode, param, f1, f2]

    for k = 1:size(paretoPoints,1)
        p1 = paretoPoints(k,1);
        p2 = paretoPoints(k,2);
        f1 = paretoPoints(k,3);
        f2 = paretoPoints(k,4);
        color = [0 0.4470 0.7410];
        scatter(f1, f2, 60, color, 'filled');
        text(f1+0.001, f2, sprintf('%.2f; %.2f; %.2f\n', p1, p2, 1-p1-p2), 'FontSize',8, 'Color', color);
        fprintf('%.2f; %.2f; %.2f\n', p1, p2, 1-p1-p2)
    end
    hold off

    subplot(3,1,3)
    hold on
    title("gamma")
    grid on;
    obj = [];
    c = 1;
    for p3_last = linspace(0,1,11)
        for i = 1:numel(params_used(:,1))
            p3 = 1 - params_used(i,1) - params_used(i,2);
            if p3 == p3_last
                obj = [obj; obj1(i) obj2(i)];
            end
        end
        color = colors(c,:);
        if ~isempty(obj)
            plot(obj(:,1), obj(:,2), "o-",'Color',color,'MarkerSize',4, 'MarkerFaceColor', color)
        end
        obj = [];
        c = min(c+1, 12);
    end

    % pareto points big dark blue
    paretoPoints = results(is_pareto, :);  % [mode, param, f1, f2]

    for k = 1:size(paretoPoints,1)
        p1 = paretoPoints(k,1);
        p2 = paretoPoints(k,2);
        f1 = paretoPoints(k,3);
        f2 = paretoPoints(k,4);
        color = [0 0.4470 0.7410];
        scatter(f1, f2, 60, color, 'filled');
        text(f1+0.001, f2, sprintf('%.2f; %.2f; %.2f\n', p1, p2, 1-p1-p2), 'FontSize',8, 'Color', color);
        fprintf('%.2f; %.2f; %.2f\n', p1, p2, 1-p1-p2)
    end
    hold off
    
    xlabel('Cost');
    ylabel('Trivial Sol.');
    % legend('All Samples','Pareto','Location','best');
    sgtitle('Pareto Front for Scale Free');
    savefig(gcf, "figures_new/pareto_ScaleFree_" + ddp)
    print(gcf, "figures_new/pareto_ScaleFree_" + ddp + ".eps", "-depsc")
    
    % pareto front (alpha; beta)
    % 0.10; 0.00
    % 0.90; 0.00
    % 1.00; 0.00
    % more exact fixed beta = 0;
    % 0.01; 0.00
end

function [cost, trivial] = DDP_sfg(alpha, beta, ddp, n)
    disp("alpha: " + alpha + ", beta: " + beta)
    gamma = 1 - alpha - beta;
    % disp("alpha = " + alpha + ", beta = " + beta + ", gamma = " + gamma)
    fract_targ = 0.1;
    fract_dist = 0.1;
    params = num2cell([n, alpha, beta, gamma, 1, 1]);
    [costs, ~, trivials] = run_test( ...
        @sfg, ...
        params, ...
        "sample_size", 500, ...
        "fraction_targets", fract_targ, ...
        "fraction_disturbances", fract_dist, ...
        "ddp", ddp, ...
        "seed", 1000 ...
        );
    cost = mean(costs);
    trivial = mean(trivials);
end

