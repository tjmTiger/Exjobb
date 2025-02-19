clear all;
clear;
clc;
% G = ER_Graph(100, 0.1, 0, 100);
% G = WattsStrogatz(100, 4, 0);
% G = SFG_dir(100, 0.3, 0.5, 0.2, 0);

% figure(1);
% plot(G)
% plot(G,'NodeColor','k','Layout','circle');

load 'konect.mat' data;
for i = 1:27
    G = cell2mat(data(i));
    figure(i);
    spy(G)
end