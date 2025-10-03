clear;clc;close all;
% figure;
% b = [];
% k_r = [];
% for beta = [0 0.44 0.62 0.72 0.78 0.83]
%     n = 50;
%     alpha = (1 - beta)/2; gamma = alpha;
%     delta_in = 1; delta_out = 1;
%     k_real = [];
%     disp(alpha + ", " + beta + ", " + gamma)
%     for seed = 1:200
%         G = sfg(n, alpha, beta, gamma, delta_in, delta_out, seed);
%         k_real(end+1) = 2*numedges(G) / numnodes(G);
%     end
%     b(end+1) = beta;
%     k_r(end+1) = mean(k_real);
% end
% hold on
% plot(b, k_r, "-o")
% hold off
% grid on
% xlabel("beta")
% ylabel("k_{real}")
% 
% %%
% k_e = [];
% k_r = [];
% for k_expected = 1:1:7
%     n = 100;
%     p = (k_expected)/(2*(n-1));
%     disp(p)
%     k_real = [];
%     for seed = 1:500
%         G = erdos_renyi(n, p, seed);
%         k_real(end+1) = 2*numedges(G) / numnodes(G);
%     end
%     k_e(end+1) = k_expected;
%     k_r(end+1) = mean(k_real);
% end
% figure;
% hold on
% plot(k_e, k_r, "-o")
% hold off
% xlabel("k_{expected}")
% ylabel("k_{real}")
% %%
% clear;clc;
% k_expected = 4;
% n = 10;
% k_half = floor(k_expected/2);
% p_ws = 0;
% k_real = [];
% % G = watts_strogatz(n, k_half, p_ws,1);
% for seed = 1:500
%     G = watts_strogatz(n, k_half, p_ws,seed);
%     k_real(end+1) = 2*numedges(G) / numnodes(G);
% end
% figure;plot(G)
% k_r = mean(k_real);
% disp("k_e: " + k_expected + ", k_r: " + k_r)
% %%
% G = [0 1; 0 0];
% G = digraph(G);
% figure;plot(G)
% disp(2*numedges(G) / numnodes(G))
% %%
% clear;clc;close all;
% n = 11;
% beta = 0;
% alpha = (1 - beta)/2;
% gamma = alpha;
% delta_in = 1; delta_out = 1;
% G = sfg(n, alpha, beta, gamma, delta_in, delta_out, 1);
% figure;plot(G)
% disp(2*numedges(G) / numnodes(G))