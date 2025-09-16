clear; clc; close all;
% multithreading()

%-----------------------------------------------%
%                                               %
%                 Optimal for SF                %
%                                               %
%-----------------------------------------------%
% min cost
% [0.03, 2, 0.9, 0.02, 0, 0.98] p, k, pw, alpha, beta, gamma

% % state feedback
% TEST(@test_fraction_size, "Size", 100:50:200, "Fractions", 0.05:0.05:0.3, "sample_size", 500)
% TEST(@test_different_fractions, "FracDist", 0.05:0.05:0.3, "FracTarg", 0.05:0.05:0.3, "sample_size", 500)
% TEST(@test_size, "", "", "Size", 30:30:200, "fract_targ", 0.3, "sample_size", 500, "boxchart", true)
% TEST(@test_node_degree, "", "", "AverageDegree", 2:1:7, "sample_size", 500, "boxchart", true)
% TEST(@test_strogatz_fraction_beta, "Fractions", 0.05:0.05:0.3, "Beta", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500)
% TEST(@test_strogatz_size_beta, "Size", 100:50:200, "Beta", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500)
% TEST(@test_scale_fraction_alphagamma, "Fractions", 0.05:0.05:0.3, "Alpha", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
% TEST(@test_scale_fraction_beta, "Fractions", 0.05:0.05:0.3, "Beta", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
% TEST(@test_scale_size_alphagamma, "Size", 100:50:200, "Alpha", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
% TEST(@test_scale_size_beta, "Size", 100:50:200, "Beta", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500)

% % output feedback
% TEST(@test_fraction_size, "Size", 100:50:200, "Fractions", 0.05:0.05:0.3, "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_different_fractions, "FracDist", 0.05:0.05:0.3, "FracTarg", 0.05:0.05:0.3, "splotample_size", 500, "ddp", "output_feedback")
% TEST(@test_size, "", "", "Size", 30:30:200, "fract_targ", 0.3, "sample_size", 500, "ddp", "output_feedback", "boxchart", true)
% TEST(@test_node_degree, "", "", "AverageDegree", 2:1:7, "sample_size", 500, "ddp", "output_feedback", "boxchart", true)
% TEST(@test_strogatz_fraction_beta, "Fractions", 0.05:0.05:0.3, "Beta", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_strogatz_size_beta, "Size", 100:50:200, "Beta", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_fraction_alphagamma, "Fractions", 0.05:0.05:0.3, "Alpha", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_fraction_beta, "Fractions", 0.05:0.05:0.3, "Beta", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_size_alphagamma, "Size", 100:50:200, "Alpha", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_size_beta, "Size", 100:50:200, "Beta", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")

%-----------------------------------------------%
%                                               %
%                 Optimal for OF                %
%                                               %
%-----------------------------------------------%
% min cost
% [0.03, 2, 0.8, 0.01, 0, 0.99] p, k, pw, alpha, beta, gamma

% TEST(@test_OFDF_size, "Feedback", ["output_feedback", "dynamical_feedback"], "Size", 20:10:80, "fract_targ", 0.3, "sample_size", 500)
% TEST(@test_OFDF_node_degree, "Feedback", ["output_feedback", "dynamical_feedback"], "AverageDegree", 2:1:7, "sample_size", 500, "boxchart", true, "size", 50)
% TEST(@test_OFDF_strogatz_beta, "Feedback", ["output_feedback", "dynamical_feedback"], "Beta", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500, "size", 50)
% TEST(@test_OFDF_scale_alphagamma, "Feedback", ["output_feedback", "dynamical_feedback"], "Alpha", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "size", 50)
% TEST(@test_OFDF_scale_beta, "Feedback", ["output_feedback", "dynamical_feedback"], "Beta", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "size", 50)
