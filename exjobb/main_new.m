clear; clc;


sample_size = 1; % sample size of graph
%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

% fraction and size
t = Tests();
for graph_generating_algorithm = {@erdos_renyi, @watts_strogatz, @sfg}
    for n = 100:50:200
        for fraction = 0.05:0.05:0.3
            % test = Test(@erdos_renyi, [n, 2, 0.2], "number_of_tests", sample_size);
            t.run(@erdos_renyi, [n, 2, 0.2], "number_of_tests", sample_size)
            disp("done")
        end
    end
end