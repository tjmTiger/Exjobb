clear; clc;

% Use Parallel pools to speed up the code.
p = gcp("nocreate");
if isempty(p) % If no pool, create new one.
    parpool(6)
else % If there is a pool, delete it and create new one.
    delete(gcp('nocreate'))
    parpool(6)
end

%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

% fraction and size
tic
t = Tests();
for graph_generating_algorithm = {@erdos_renyi, @watts_strogatz, @sfg}
    for n = 100:50:200
        for fraction = 0.05:0.05:0.3
            disp(graph_generating_algorithm)
            t.run(@erdos_renyi, [n, 2, 0.2], "sample_size", 200, "fraction_targets", 0.1)
        end
    end
end
t.plot()
toc