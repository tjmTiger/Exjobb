clear;
close all;
clc;


fig = openfig("figures/erdos_fract.fig");
ax_handles = findall(fig, 'Type', 'Axes');
% ax_handles(1).YLim = [-1, 1];
ax_handles(1).FontSize = 10;
ax_handles(2).FontSize = 10;
ax_handles(3).FontSize = 10;

% figure();
% s1 = subplot(3,1,1);
% s2 = subplot(3,1,2);
% s3 = subplot(3,1,3);
% 
% fig1 = get(ax_handles(1),'children');
% copyobj(fig1, s1)


