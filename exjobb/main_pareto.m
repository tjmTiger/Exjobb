clear; clc;
multithreading()
% todo: pareto optimization => use gamultiobj()

f = @(x)DDP_erods_renyi(x(1), x(2), x(3));

x = gamultiobj(f, 3, [0,1,1], 0.9, [], [], 0, 1);

plot(x(:,1),x(:,2),'ko')
xlabel('x(1)')
ylabel('x(2)')
title('Pareto Points in Parameter Space')


function [cost, runtime] = DDP_erods_renyi(p, fract_targ, fract_dist)
    n = 100;
    G = erdos_renyi(n, p, randi(500));
    params = num2cell([n, p]);
    [costs, runtimes, ~] = run_test( ...
        @erdos_renyi, ...
        params, ...
        "sample_size", 1, ...
        "fraction_targets", fract_targ, ...
        "fraction_disturbances", fract_dist, ...
        "ddp", "state_feedback" ...
        );
    cost = mean(costs);
    runtime = mean(runtimes);
end
