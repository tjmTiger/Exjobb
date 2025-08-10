clear; clc;
multithreading()

test_fraction_size("sample_size", 200)
test_different_fractions("sample_size", 200)
test_size("sample_size", 200)
test_node_degree("sample_size", 200)
test_strogatz_fraction_beta("sample_size", 200)

% To do:
% Watts Strogatz
%   - fractions beta
%   - size beta
% Scale Free
%   - fractions alpha-gamma
%   - fractions beta
%   - size alpha-gamma
%   - size beta