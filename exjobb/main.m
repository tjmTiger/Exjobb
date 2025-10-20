clear; clc; close all;
multithreading()

%-----------------------------------------------%
%                                               %
%                 Optimal for SF                %
%                                               %
%-----------------------------------------------%
% % state feedback
% TEST(@test_fraction_size, "Size", 100:50:200, "Fractions", 0.05:0.05:0.3, "sample_size", 500)
% TEST(@test_different_fractions, "FracDist", 0.05:0.05:0.3, "FracTarg", 0.05:0.05:0.3, "sample_size", 500)
% TEST(@test_size, "", "", "Size", 30:30:200, "fract_targ", 0.3, "sample_size", 500, "boxchart", true)
% TEST(@test_node_degree, "", "", "AverageDegree", 2:1:7, "sample_size", 500, "boxchart", true)
% TEST(@test_erdos_fraction_p, "Fractions", 0.05:0.05:0.3, "$p$", [0.03 0.2 0.4 0.6 0.8 1], "graph_generating_algorithm", "Erdos Renyi", "sample_size", 500)
% TEST(@test_strogatz_fraction_pws, "Fractions", 0.05:0.05:0.3, "$p_{ws}$", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500)
% TEST(@test_strogatz_size_pws, "Size", 100:50:200, "$p_{ws}$", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500)
% TEST(@test_scale_fraction_alphagamma, "Fractions", 0.05:0.05:0.3, "$\alpha$", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
% TEST(@test_scale_fraction_beta, "Fractions", 0.05:0.05:0.3, "$\beta$", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
% TEST(@test_scale_size_alphagamma, "Size", 100:50:200, "$\alpha$", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
% TEST(@test_scale_size_beta, "Size", 100:50:200, "$\beta$", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500)

% % output feedback
% TEST(@test_fraction_size, "Size", 100:50:200, "Fractions", 0.05:0.05:0.3, "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_different_fractions, "FracDist", 0.05:0.05:0.3, "FracTarg", 0.05:0.05:0.3, "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_size, "", "", "Size", 30:30:200, "fract_targ", 0.3, "sample_size", 500, "ddp", "output_feedback", "boxchart", true)
% TEST(@test_node_degree, "", "", "AverageDegree", 2:1:7, "sample_size", 500, "ddp", "output_feedback", "boxchart", true)
% TEST(@test_erdos_fraction_p, "Fractions", 0.05:0.05:0.3, "$p$", [0.03 0.2 0.4 0.6 0.8 1], "graph_generating_algorithm", "Erdos Renyi", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_strogatz_fraction_pws, "Fractions", 0.05:0.05:0.3, "$p_{ws}$", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_strogatz_size_pws, "Size", 100:50:200, "$p_{ws}$", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_fraction_alphagamma, "Fractions", 0.05:0.05:0.3, "$\alpha$", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_fraction_beta, "Fractions", 0.05:0.05:0.3, "$\beta$", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_size_alphagamma, "Size", 100:50:200, "$\alpha$", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")
% TEST(@test_scale_size_beta, "Size", 100:50:200, "$\beta$", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "ddp", "output_feedback")

%-----------------------------------------------%
%                                               %
%                 Optimal for OF                %
%                                               %
%-----------------------------------------------%
% default_p = [0.1, 2, 0.9, 1, 0, 0]; % p, k, pw, alpha, beta, gamma

% TEST(@test_OFDF_node_degree, "Feedback", ["output_feedback", "dynamical_feedback"], "AverageDegree", 2:1:7, "sample_size", 500, "size", 50, "default_p", default_p, "ddp", "dynamical_feedback")
% TEST(@test_OFDF_erdos_p, "Feedback", ["output_feedback", "dynamical_feedback"], "$p$", [0.1 0.2 0.4 0.6 0.8 1], "graph_generating_algorithm", "Erdos Renyi", "sample_size", 500, "size", 50, "default_p", default_p, "ddp", "dynamical_feedback")
% TEST(@test_OFDF_strogatz_pws, "Feedback", ["output_feedback", "dynamical_feedback"], "$p_{ws}$", 0:0.1:1, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500, "size", 50, "default_p", default_p, "ddp", "dynamical_feedback")
% TEST(@test_OFDF_scale_alphagamma, "Feedback", ["output_feedback", "dynamical_feedback"], "$\alpha$", 0:0.1:1, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "size", 50, "default_p", default_p, "ddp", "dynamical_feedback")
% TEST(@test_OFDF_scale_beta, "Feedback", ["output_feedback", "dynamical_feedback"], "$\beta$", 0:0.1:0.8, "graph_generating_algorithm", "Scale Free", "sample_size", 500, "size", 50, "default_p", default_p, "ddp", "dynamical_feedback")
% TEST(@test_OFDF_size, "Feedback", ["output_feedback", "dynamical_feedback"], "Size", 20:10:50, "fract_targ", 0.3, "sample_size", 500, "default_p", default_p, "ddp", "dynamical_feedback")
