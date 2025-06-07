clear;
close all;
clc;

mycolors = [
238,64,53;
243,155,54;
123,192,67;
3,146,207;
17,0,255;
175,56,255;
]./255;

n_graphs = 200; % sample size of graph
%-----------------------------------------------%
%                                               %
%     Distrubance and target node fractions     %
%                                               %
%-----------------------------------------------%

for graph_name = ["Erdos Renyi", "Watts Strogatz", "Scale Free"]
    figure();
    for n = 100:50:200
        display_name = string(n);
        results_all = [];
        results_time_all = [];
        results_trivial_all = [];
        fract_T_D = 0.05:0.05:0.3;
        for j = fract_T_D
            fract_targ = j;
            fract_dist = j;
            results = zeros(1, n_graphs);
            results_time = zeros(1, n_graphs);
            results_trivial = zeros(1, n_graphs);
            for i = 1:n_graphs
                if graph_name == "Erdos Renyi"
                    p = 1.5*(log10(n)/n);
                    G = ER_Graph(n, p, 1, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("Erdos: " + i)
                elseif graph_name == "Watts Strogatz"
                    k = 2;
                    beta = 0.2;
                    G = WattsStrogatz(n, k, beta, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("Strogatz: " + i)
                else
                    alpha = 0.4;
                    beta = 0.2;
                    gamma = 0.4;
                    G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("SFG: " + i)
                end
            end
            results_all(end+1,:) = results;
            results_time_all(end+1,:) = results_time;
            results_trivial_all(end+1,:) = results_trivial;
        end
        subplot(1,3,1);
        hold on;
        
        mean_list = [];
        for r = 1:size(results_all, 1)
            mean_list(:,end+1) = mean(results_all(r,:));
        end
        plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
        title("Cost")
        ylabel("Cost [-]")
        xlabel("Fractions")
        xlim([min(fract_T_D) max(fract_T_D)])
        ax = gca; 
        ax.ColorOrder = mycolors;
        
        hold off;
    
        subplot(1,3,2);
        hold on;
        
        mean_list = [];
        for r = 1:size(results_time_all, 1)
            mean_list(:,end+1) = mean(results_time_all(r,:));
        end
        plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
        title(graph_name + newline + "Runtime")
        ylabel("Time [s]")
        xlabel("Fractions")
        xlim([min(fract_T_D) max(fract_T_D)])
        ax = gca; 
        ax.ColorOrder = mycolors;
    
        hold off;
    
        subplot(1,3,3);
        hold on;
        
        mean_list = [];
        for r = 1:size(results_trivial_all, 1)
            mean_list(:,end+1) = mean(results_trivial_all(r,:)); 
        end
        plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
        title("Trivial solutions")
        ylabel("Index [-]")
        xlabel("Fractions")
        ylim([0 1])
        xlim([min(fract_T_D) max(fract_T_D)])
        ax = gca; 
        ax.ColorOrder = mycolors;
        lgd = legend;
        title(lgd,'size')
        hold off;
    end
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    saveas(gcf, "figures_new/" + graph_name + " size_fract.jpg")
end
%
%-----------------------------------------------%
%                                               %
%               Different Fractions             %
%                                               %
%-----------------------------------------------%

for graph_name = ["Erdos Renyi", "Watts Strogatz", "Scale Free"]
    figure();
    n = 100;
    fract_T_D = 0.05:0.05:0.3;
    for fract_dist = fract_T_D
        display_name = string(fract_dist);
        results_all = [];
        results_time_all = [];
        results_trivial_all = [];
        for fract_targ = fract_T_D
            results = zeros(1, n_graphs);
            results_time = zeros(1, n_graphs);
            results_trivial = zeros(1, n_graphs);
            for i = 1:n_graphs
                if graph_name == "Erdos Renyi"
                    p = 1.5*(log10(n)/n);
                    G = ER_Graph(n, p, 1, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("Erdos: " + i)
                elseif graph_name == "Watts Strogatz"
                    k = 2;
                    beta = 0.2;
                    G = WattsStrogatz(n, k, beta, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("Strogatz: " + i)
                else
                    alpha = 0.4;
                    beta = 0.2;
                    gamma = 0.4;
                    G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("SFG: " + i)
                end
            end
            results_all(end+1,:) = results;
            results_time_all(end+1,:) = results_time;
            results_trivial_all(end+1,:) = results_trivial;
        end
        subplot(1,3,1);
        hold on;
        
        mean_list = [];
        for r = 1:size(results_all, 1)
            mean_list(:,end+1) = mean(results_all(r,:));
        end
        plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
        title("Cost")
        ylabel("Cost [-]")
        xlabel("Targ.Frac.")
        xlim([min(fract_T_D) max(fract_T_D)])
        ax = gca; 
        ax.ColorOrder = mycolors;
        
        hold off;
    
        subplot(1,3,2);
        hold on;
        
        mean_list = [];
        for r = 1:size(results_time_all, 1)
            mean_list(:,end+1) = mean(results_time_all(r,:));
        end
        plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
        title(graph_name + newline + "Runtime")
        ylabel("Time [s]")
        xlabel("Targ.Frac.")
        xlim([min(fract_T_D) max(fract_T_D)])
        ax = gca; 
        ax.ColorOrder = mycolors;
    
        hold off;
    
        subplot(1,3,3);
        hold on;
        
        mean_list = [];
        for r = 1:size(results_trivial_all, 1)
            mean_list(:,end+1) = mean(results_trivial_all(r,:)); 
        end
        plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
        title("Trivial solutions")
        ylabel("Index [-]")
        xlabel("Targ.Frac.")
        ylim([0 1])
        xlim([min(fract_T_D) max(fract_T_D)])
        ax = gca; 
        ax.ColorOrder = mycolors;
        lgd = legend;
        title(lgd,'Dist.Frac.')
        hold off;
    end
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    saveas(gcf, "figures_new/" + graph_name + " diff_fract.jpg")
end

%%
%-----------------------------------------------%
%                                               %
%                   Graph Size                  %
%                                               %
%-----------------------------------------------%

fract_targ = 0.1;
fract_dist = 0.1;

for graph_name = ["Erdos Renyi", "Watts Strogatz", "Scale Free"]
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    graph_size = 30:30:180;
    for j = graph_size
        n = j;
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs
            if graph_name == "Erdos Renyi"
                p = 0.03; % 1.5*(log10(n)/n);
                G = ER_Graph(n, p, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Erdos: " + i)
            elseif graph_name == "Watts Strogatz"
                k = 2;
                beta = 0.2;
                G = WattsStrogatz(n, k, beta, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Strogatz: " + i)
            else
                alpha = 0.4;
                beta = 0.2;
                gamma = 0.4;
                G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("SFG: " + i)
            end
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end

    figure();
    subplot(1,3,1);
    hold on;

    for r = 1:size(results_all, 1)
        add2boxchart(results_all(r,:), string(graph_size(r)), "Cost", "Cost [-]", "Size")
    end

    hold off;

    subplot(1,3,2);
    hold on;

    for r = 1:size(results_time_all, 1)
        add2boxchart(results_time_all(r,:), string(graph_size(r)), graph_name + newline + "Runtime", "Time [s]", "Size")
    end

    hold off;

    subplot(1,3,3);
    hold on;

    for r = 1:size(results_trivial_all, 1)
        add2boxchart(results_trivial_all(r,:), string(graph_size(r)), "Trivial solutions", "Index [-]", "Size")
    end
    ylim([0 1])
    hold off;

    saveas(gcf, "figures_new/" + graph_name + " size.jpg")
end
%%
%-----------------------------------------------%
%                                               %
%                 Connectivity                  %
%                                               %
%-----------------------------------------------%

fract_targ = 0.3;
fract_dist = 0.1;
n = 100;
error = 0;
for graph_name = ["Erdos Renyi", "Watts Strogatz"]
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    node_degree = 2:1:7;
    for j = node_degree
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs
            if graph_name == "Erdos Renyi"
                p = 2.2*j/(n-1); % 1.5*(log10(n)/n);
                G = ER_Graph(n, p, 1, i);
                error = error + (numedges(G)/numnodes(G)-j);
                [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                disp("Erdos: " + i)
            elseif graph_name == "Watts Strogatz"
                if mod(j,2) == 0
                    k = floor(j/2);
                    beta = 0.2;
                    G = WattsStrogatz(n, k, beta, i);
                    [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
                    disp("Strogatz: " + i)
                else
                    results(i) = NaN;
                    results_time(i) = NaN;
                    results_trivial(i) = NaN;
                end
            end
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end

    figure();
    subplot(1,3,1);
    hold on;

    for r = 1:size(results_all, 1)
        add2boxchart(results_all(r,:), string(node_degree(r)), "Cost", "Cost [-]", "Average degree")
    end

    hold off;

    subplot(1,3,2);
    hold on;

    for r = 1:size(results_time_all, 1)
        add2boxchart(results_time_all(r,:), string(node_degree(r)), graph_name + newline + "Runtime", "Time [s]", "Average degree")
    end

    hold off;

    subplot(1,3,3);
    hold on;

    for r = 1:size(results_trivial_all, 1)
        add2boxchart(results_trivial_all(r,:), string(node_degree(r)), "Trivial solutions", "Index [-]", "Average degree")
    end
    ylim([0 1])
    hold off;

    saveas(gcf, "figures_new/" + graph_name + " node_degree.jpg")
end

disp("Mean node degree error in Erdos Renyi: " + error/(n_graphs*numel(node_degree))) % -0.0096

%% scale free
tot_mean_degree = [];
results = zeros(1, n_graphs);
todo = n_graphs*6;
k = 2; % start value for node degree
figure();
for beta = [0.55 0.70 0.79 0.845 0.871 0.893]
    for i = 1:n_graphs
        % beta = 0.1;
        alpha = (1-beta)/2;
        gamma = (1-beta)/2;
        delta_in = 10;
        delta_out = 10;
        G = SFG_dir(n, alpha, beta, gamma, delta_in, delta_out, i);
        [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist); 
        todo = todo-1;
        disp("Left: " + todo)
    end
    subplot(1,3,1);
    hold on
    add2boxchart(results, string(k), "Cost", "Cost [-]", "average degree")
    hold off
    subplot(1,3,2);
    hold on
    add2boxchart(results_time, string(k), "Scale free" + newline + "Runtime", "Time [s]", "average degree")
    hold off
    subplot(1,3,3);
    hold on
    add2boxchart(results_trivial, string(k), "Trivial solutions", "Index [-]", "average degree")
    hold off
    ylim([0 1])
    k = k+1;
end
saveas(gcf, "figures_new/" + "Scale Free" + " node_degree.jpg")

%%
%-----------------------------------------------%
%                                               %
%                  Extra tests                  %
%                   Scale Free                  %
%                                               %
%-----------------------------------------------%
%% fractions alpha - gamma

graph_name = "Scale Free";
figure();
n = 100;
for alpha = 0.2:0.2:0.8
    display_name = string(alpha);
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    fract_T_D = 0.05:0.05:0.3;
    for j = fract_T_D
        fract_targ = j;
        fract_dist = j;
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs     
            beta = 0.2;
            gamma = 0.8-alpha;
            G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
            [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
            disp("SFG: " + i)
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end
    subplot(1,3,1);
    hold on;
    
    mean_list = [];
    for r = 1:size(results_all, 1)
        mean_list(:,end+1) = mean(results_all(r,:));
    end
    plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
    title("Cost")
    ylabel("Cost [-]")
    xlabel("Fractions")
    xlim([min(fract_T_D) max(fract_T_D)])
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    hold off;

    subplot(1,3,2);
    hold on;
    
    mean_list = [];
    for r = 1:size(results_time_all, 1)
        mean_list(:,end+1) = mean(results_time_all(r,:));
    end
    plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
    title(graph_name + newline + "Runtime")
    ylabel("Time [s]")
    xlabel("Fractions")
    xlim([min(fract_T_D) max(fract_T_D)])
    ax = gca; 
    ax.ColorOrder = mycolors;

    hold off;

    subplot(1,3,3);
    hold on;
    
    mean_list = [];
    for r = 1:size(results_trivial_all, 1)
        mean_list(:,end+1) = mean(results_trivial_all(r,:)); 
    end
    plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
    title("Trivial solutions")
    ylabel("Index [-]")
    xlabel("Fractions")
    ylim([0 1])
    xlim([min(fract_T_D) max(fract_T_D)])
    ax = gca; 
    ax.ColorOrder = mycolors;
    lgd = legend;
    title(lgd,'alpha')
    hold off;
end
fontsize(12,"points")
position = get(gcf, 'Position');
position = [100, 100, 600, 600];
saveas(gcf, "figures_new/" + graph_name + " alpha_fract.jpg")

%% fractions beta
graph_name = "Scale Free";
figure();
n = 100;
for beta = 0.2:0.2:0.8
    display_name = string(beta);
    results_all = [];
    results_time_all = [];
    results_trivial_all = [];
    fract_T_D = 0.05:0.05:0.3;
    for j = fract_T_D
        fract_targ = j;
        fract_dist = j;
        results = zeros(1, n_graphs);
        results_time = zeros(1, n_graphs);
        results_trivial = zeros(1, n_graphs);
        for i = 1:n_graphs     
            alpha = (1-beta)/2;
            gamma = 1-beta-alpha;
            G = SFG_dir(n, alpha, beta, gamma, 1, 1, i);
            [results(i), results_time(i), results_trivial(i)] = decouple(G, fract_targ, fract_dist);
            disp("SFG: " + i)
        end
        results_all(end+1,:) = results;
        results_time_all(end+1,:) = results_time;
        results_trivial_all(end+1,:) = results_trivial;
    end
    subplot(1,3,1);
    hold on;
    
    mean_list = [];
    for r = 1:size(results_all, 1)
        mean_list(:,end+1) = mean(results_all(r,:));
    end
    plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
    title("Cost")
    ylabel("Cost [-]")
    xlabel("Fractions")
    xlim([min(fract_T_D) max(fract_T_D)])
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    hold off;

    subplot(1,3,2);
    hold on;
    
    mean_list = [];
    for r = 1:size(results_time_all, 1)
        mean_list(:,end+1) = mean(results_time_all(r,:));
    end
    plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
    title(graph_name + newline + "Runtime")
    ylabel("Time [s]")
    xlabel("Fractions")
    xlim([min(fract_T_D) max(fract_T_D)])
    ax = gca; 
    ax.ColorOrder = mycolors;

    hold off;

    subplot(1,3,3);
    hold on;
    
    mean_list = [];
    for r = 1:size(results_trivial_all, 1)
        mean_list(:,end+1) = mean(results_trivial_all(r,:)); 
    end
    plot(fract_T_D, mean_list, "-o", 'DisplayName', display_name)
    title("Trivial solutions")
    ylabel("Index [-]")
    xlabel("Fractions")
    ylim([0 1])
    xlim([min(fract_T_D) max(fract_T_D)])
    ax = gca; 
    ax.ColorOrder = mycolors;
    lgd = legend;
    title(lgd,'beta')
    hold off;
end
fontsize(12,"points")
position = get(gcf, 'Position');
position = [100, 100, 600, 600];
saveas(gcf, "figures_new/" + graph_name + " beta_fract.jpg")

%-----------------------------------------------%
%                                               %
%                   Functions                   %
%                                               %
%-----------------------------------------------%

function add2boxchart(results, test_name, title_name, ylabel_name, xlabel_name)
    fontsize(12,"points")
    position = get(gcf, 'Position');
    position = [100, 100, 600, 600];
    boxchart(categorical(1:numel(results), 1:numel(results), repmat(test_name, 1, numel(results))), results)
    if nargin >= 3
        title(title_name)
    end
    if nargin == 4
        ylabel(ylabel_name)
        xlabel('Tests')
    elseif nargin == 5
        ylabel(ylabel_name)
        xlabel(xlabel_name)
    else
        ylabel('Cost [-]')
        xlabel('Tests')
    end
end