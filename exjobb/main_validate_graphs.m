clear; clc;

graph_sizes = 100:100:400;
sample_size = 200;

for graph_gen_alg = ["erdos renyi"] % ["erdos renyi" , "watts strogatz ring", "watts strogatz", "scalefree1", "scalefree2"]
    results = [];
    subplot_num = 1;
    figure();
    for n = graph_sizes
        degrees = [];% zeros(1,n_real);
        subplot(2,2,subplot_num)
        subplot_num = subplot_num + 1;
        hold on
        title(graph_gen_alg + " n = " + string(n))
        for sample = 1:sample_size
            switch graph_gen_alg
                case "erdos renyi"
                    p = 0.03;
                    G = erdos_renyi(n, p, sample);
                case "watts strogatz ring"
                    k = 2;
                    beta = 0;
                    G = watts_strogatz(n, k, beta, sample);
                case "watts strogatz"
                    k = 2;
                    beta = 0.2;
                    G = watts_strogatz(n, k, beta, sample);
                case "scalefree1"
                    alpha = 0.5;
                    beta = 0;
                    gamma = 0.5;
                    G = sfg(n, alpha, beta, gamma, 1, 1, sample);
                case "scalefree2"
                    alpha = 0.1;
                    beta = 0.8;
                    gamma = 0.1;
                    G = sfg(n, alpha, beta, gamma, 1, 1, sample);
            end
            for i = 1:numnodes(G)
                degrees(end+1) = indegree(G, i);%  + outdegree(G, i);
            end
        end

        switch graph_gen_alg
            case {"scalefree1", "scalefree2"}
                deg = degrees;                       % degree of each node
                [k_counts, k_bins] = histcounts(deg, 'BinMethod', 'integers');
                 
                % Shift bin centers
                k_vals = k_bins(1:end-1) + diff(k_bins)/2;
                 
                % Remove zero-count bins (for log scale)
                mask = (k_counts > 0);
                k_vals = k_vals(mask);
                k_counts = k_counts(mask);
                loglog(k_vals, k_counts./sample_size)
            otherwise
                % yyaxis left
                [counts,bins] = histcounts(degrees);
                histogram('BinEdges',bins,'BinCounts',counts/sample_size)
        end

        xlabel("Degree")
        ylabel("Frequency")
        xlim([0,max(degrees)])

        switch graph_gen_alg
            case "erdos renyi"
                k = 0:nchoosek(n, 2); % n-1;      % Possible degree values (from 0 to n-1)
                P_k = binopdf(k, n, p);  % Binomial degree distribution
                plot(k, n.*P_k, "green");

                % lambda = (n-1)*p;
                % P_poisson = poisspdf(k, lambda);
                % plot(k, n.*P_poisson, "red");
            case "watts strogatz"
                % pass
            case {"scalefree1", "scalefree2"}
                % todo: try with polyfit
                k = 0:nchoosek(n, 2);

                X_in = 1 + (1+1*(alpha+gamma)) / (alpha + beta);
                % X_out = 1 + (1+1*(alpha+gamma)) / (gamma + beta);

                P_k_in = (1/zeta(X_in)).*k.^-X_in;
                % P_k_out = (1/zeta(X_out)).*k.^-X_out;

                P_k = P_k_in;
                % yyaxis right
                loglog(k, n.*P_k, "r");
                set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')
        end
        hold off
    end
    interpreter = 'latex';
    all_text = findall(gcf, "-property", "Interpreter");
    set(all_text, "Interpreter", interpreter)

    all_text = findall(gcf, "-property", "TickLabelInterpreter");
    set(all_text, "TickLabelInterpreter", interpreter)

    all_text = findall(gcf, "-property", "Fontsize");
    set(all_text, "Fontsize", 12) % 12

    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    savefig(gcf, "figures_new/validation_" + erase(graph_gen_alg," "))
    print(gcf, "figures_new/validation_" + erase(graph_gen_alg," ") + ".eps", "-depsc")
end