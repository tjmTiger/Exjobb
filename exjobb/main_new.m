clear; clc;


sample_size = 10; % sample size of graph
%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

% fraction and size
for graph_generating_algorithm = {@erdos_renyi, @watts_strogatz, @sfg}
    for n = 100:50:200
        for fraction = 0.05:0.05:0.3
            test = Test(@WattsStrogatz, [n, 2, 0.2], "number_of_tests", sample_size);
            test.run()
        end
    end
end