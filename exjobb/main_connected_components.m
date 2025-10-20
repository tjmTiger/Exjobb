clear;clc;
sample_size = 500;

figure;
i = 1;
% each size in own subplot
for n = [20 50 100 180]
    subplot(2,2,i); i = i+1;
    title("Watts Strogatz $n$ = " + n)
    % each k_ws new entry in legend
    for k_ws = 2:2:4
        points = [];
        % p_ws on x-axis
        for p_ws = 0:0.1:1
            n_comp = zeros(sample_size,1);
            for seed = 1:sample_size
                G = watts_strogatz(n,k_ws/2,p_ws,seed);
                n_comp(seed) = get_connected_components(G);
            end
            points(end+1) = mean(n_comp);
        end
        hold on
        plot(0:0.1:1, points, 'DisplayName', "$k_{ws}$ = " + k_ws)
        hold off
    end
    xlabel("$p_{ws}$")
    ylabel("Components")
end

% figure;
% i = 1;
% title("Erdos Renyi")
% % each k_ws new entry in legend
% for n = [20 50 100 180]
%     points = [];
%     % p on x-axis
%     for p = 0.03:0.01:0.2
%         n_comp = zeros(sample_size,1);
%         for seed = 1:sample_size
%             G = erdos_renyi_unconnected(n,p,seed);
%             n_comp(seed) = get_connected_components(G);
%         end
%         points(end+1) = mean(n_comp);
%     end
%     hold on
%     plot(0.03:0.01:0.2, points, 'DisplayName', "$n$ = " + n)
%     hold off
% end
% xlabel("$p$")
% ylabel("Components")

% post process figure
legend;

interpreter = 'latex';
all_text = findall(gcf, "-property", "Interpreter");
set(all_text, "Interpreter", interpreter)

all_text = findall(gcf, "-property", "TickLabelInterpreter");
set(all_text, "TickLabelInterpreter", interpreter)

all_text = findall(gcf, "-property", "Fontsize");
set(all_text, "Fontsize", 12)


function n_comp = get_connected_components(G)
    Ag = full(adjacency(G));
    % Get number of components (conncomp only works on undir. graphs)
    [~, binsize] = conncomp(graph(sparse(Ag + Ag')));
    n_comp = length(binsize);
end