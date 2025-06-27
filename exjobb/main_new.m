clear; clc;

% Use Parallel pools to speed up the code.
p = gcp("nocreate");
if isempty(p) % If no pool, create new one.
    parpool(6)
else % If there is a pool, but its small, delete it and create new one.
    if p.NumWorkers < 6
        delete(gcp('nocreate'))
        parpool(6)
    end
end

%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

% fraction and size
for graph_generating_algorithm = {{@erdos_renyi, [100, 0.03]}, {@watts_strogatz, [100, 2, 0.2]}, {@sfg, [100, 0.4, 0.2, 0.4, 1, 1]}}
    t = Tests();
    figure();
    for n = 100:50:200
        for fraction = 0.05:0.05:0.3
            t.run(graph_generating_algorithm{1}{1}, graph_generating_algorithm{1}{2}, "sample_size", 10, "fraction_targets", 0.1);
        end
        t.plot(0.05:0.05:0.3)
    end
end
