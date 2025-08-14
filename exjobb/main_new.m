clear; clc;
multithreading()

% state feedback
TEST(@test_fraction_size, "Size", 100:50:200, "Fractions", 0.05:0.05:0.3, "sample_size", 500)
TEST(@test_different_fractions, "FracDist", 0.05:0.05:0.3, "FracTarg", 0.05:0.05:0.3, "sample_size", 500)
TEST(@test_size, "", "", "Size", 30:30:180, "fract_targ", 0.3, "sample_size", 500)
TEST(@test_node_degree, "", "", "AverageDegree", 2:1:7, "sample_size", 500)
TEST(@test_strogatz_fraction_beta, "Beta", 0:0.2:1, "Fractions", 0.05:0.05:0.3, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500)
TEST(@test_strogatz_size_beta, "Beta", 0:0.2:1, "Size", 30:30:180, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 500)
TEST(@test_scale_fraction_alphagamma, "Alpha", 0.2:0.2:0.8, "Fractions", 0.05:0.05:0.3, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
TEST(@test_scale_fraction_beta, "Beta", 0.2:0.2:0.8, "Fractions", 0.05:0.05:0.3, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
TEST(@test_scale_size_alphagamma, "Alpha", 0.2:0.2:0.8, "Size", 30:30:180, "graph_generating_algorithm", "Scale Free", "sample_size", 500)
TEST(@test_scale_size_beta, "Beta", 0.2:0.2:0.8, "Size", 30:30:180, "graph_generating_algorithm", "Scale Free", "sample_size", 500)


% output feedback
TEST(@test_fraction_size, "Size", 100:50:200, "Fractions", 0.05:0.05:0.3, "sample_size", 100, "ddp", "output_feedback")
TEST(@test_different_fractions, "FracDist", 0.05:0.05:0.3, "FracTarg", 0.05:0.05:0.3, "sample_size", 100, "ddp", "output_feedback")
TEST(@test_size, "", "", "Size", 30:30:180, "fract_targ", 0.3, "sample_size", 100, "ddp", "output_feedback")
TEST(@test_node_degree, "", "", "AverageDegree", 2:1:7, "sample_size", 100, "ddp", "output_feedback")
TEST(@test_strogatz_fraction_beta, "Beta", 0:0.2:1, "Fractions", 0.05:0.05:0.3, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 100, "ddp", "output_feedback")
TEST(@test_strogatz_size_beta, "Beta", 0:0.2:1, "Size", 30:30:180, "graph_generating_algorithm", "Watts Strogratz", "sample_size", 100, "ddp", "output_feedback")
TEST(@test_scale_fraction_alphagamma, "Alpha", 0.2:0.2:0.8, "Fractions", 0.05:0.05:0.3, "graph_generating_algorithm", "Scale Free", "sample_size", 100, "ddp", "output_feedback")
TEST(@test_scale_fraction_beta, "Beta", 0.2:0.2:0.8, "Fractions", 0.05:0.05:0.3, "graph_generating_algorithm", "Scale Free", "sample_size", 100, "ddp", "output_feedback")
TEST(@test_scale_size_alphagamma, "Alpha", 0.2:0.2:0.8, "Size", 30:30:180, "graph_generating_algorithm", "Scale Free", "sample_size", 100, "ddp", "output_feedback")
TEST(@test_scale_size_beta, "Beta", 0.2:0.2:0.8, "Size", 30:30:180, "graph_generating_algorithm", "Scale Free", "sample_size", 100, "ddp", "output_feedback")
