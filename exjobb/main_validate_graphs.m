clear; clc;

graph_sizes = 100:100:400;
sample_size = 200;

for graph_gen_alg = ["erdos renyi", "watts strogatz", "sfg"]
    results = [];
    subplot_num = 1;
    figure();
    for n = graph_sizes
        subplot(2,2,subplot_num)
        subplot_num = subplot_num + 1;
        hold on
        title(graph_gen_alg + " " + string(n))
        for sample = 1:sample_size
            switch graph_gen_alg
                case "erdos renyi"
                    G = erdos_renyi(n, 0.03, 1, sample);
                case "watts strogatz"
                    k = 2;
                    beta = 0.2;
                    G = watts_strogatz(n, k, beta, sample);
                case "sfg"
                    alpha = 0.4;
                    beta = 0.2;
                    gamma = 0.4;
                    G = sfg(n, alpha, beta, gamma, 1, 1, sample);
            end
            n_real = numnodes(G);
            degrees = zeros(1,n_real);
            for i = 1:n_real
                degrees(i) = indegree(G, i) + outdegree(G, i);
            end
            histogram(degrees)
        end
        hold off
    end
    
end