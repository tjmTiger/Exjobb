clear; clc;

graph_sizes = 100:100:400;
sample_size = 200;

for graph_gen_alg = ["erdos renyi", "watts strogatz ring", "watts strogatz"]%, "sfg"]
    results = [];
    subplot_num = 1;
    figure();
    for n = graph_sizes
        degrees = [];% zeros(1,n_real);
        subplot(2,2,subplot_num)
        subplot_num = subplot_num + 1;
        hold on
        title(graph_gen_alg + " " + string(n))
        for sample = 1:sample_size
            switch graph_gen_alg
                case "erdos renyi"
                    p = 0.03;
                    G = erdos_renyi(n, p, 1, sample);
                case "watts strogatz ring"
                    k = 2;
                    beta = 0;
                    G = watts_strogatz(n, k, beta, sample);
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
            for i = 1:numnodes(G)
                degrees(end+1) = indegree(G, i) + outdegree(G, i);
            end
        end
        % yyaxis left
        [counts,bins] = histcounts(degrees);
        histogram('BinEdges',bins,'BinCounts',counts/sample_size)

        xlabel("Degree")
        ylabel("Frequency")
        xlim([0,max(degrees)])

        switch graph_gen_alg
            case "erdos renyi"
                k = 0:nchoosek(n, 2); % n-1;      % Possible degree values (from 0 to n-1)

                P_k = binopdf(k, n, p);  % Binomial degree distribution
                % yyaxis right
                plot(k, n.*P_k, "green");
            
                lambda = (n-1)*p;
                P_poisson = poisspdf(k, lambda);
                % yyaxis right
                plot(k, n.*P_poisson, "red");
            case "watts strogatz"
                % pass
            case "sfg"
                k = 0:nchoosek(n, 2);

                X_in = 1 + (1+1*(alpha+gamma)) / (alpha + beta);
                X_out = 1 + (1+1*(alpha+gamma)) / (gamma + beta);
                disp(X_in + ", " + X_out)

                P_k_in = zeta(X_in).*k.^-X_in;
                % P_k_out = zeta(X_out).*k.^-X_out;

                P_k = P_k_in;
                % yyaxis right
                plot(k, n.*P_k, "r");
        end
        hold off
    end
    fontsize(10,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    print(gcf, "figures_new/validation_" + erase(graph_gen_alg," ") + ".eps", "-depsc")
end