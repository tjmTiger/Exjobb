clear; clc; close all;
multithreading()

% pareto_erdos("output_feedback", 50)
% pareto_strogatz("output_feedback", 50)
% pareto_sfg("output_feedback", 50)

% pareto_erdos("state_feedback", 100, 'params1' , linspace(0.03,0.07,17))  % (note, not worth testing p > 0.07 for SF,  ty same cost = 1 & trivial sol =  1 for all solutions)
% pareto_strogatz("state_feedback", 100)
% pareto_sfg("state_feedback", 100)