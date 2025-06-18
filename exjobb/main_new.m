clear; clc;


t = Test(@WattsStrogatz, [100, 2, 0.2], "number_of_tests", 10);

for random_graph_algorithm = [@erdos_renyi, @watts_strogatz, @sfg]
    %pass
end