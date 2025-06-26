classdef Test < handle
%TEST Stores properties and results of a test with given parameters including cost,
%runtime and ammount of trivial solutions.
%   Detailed explanation goes here

properties
    % graphs properties
    graph_algorithm 
    graph_parameters {mustBeCell}
    graph_seed {mustBeNumeric} = 0

    % test properties
    number_of_tests {mustBeNumeric} = 200
    fraction_targets {mustBeFloat} = 0.1
    fraction_disturbances {mustBeFloat} = 0.1

    % other
    results_cost
    results_time
    results_trivial

    decouple
end
methods
    function self = Test(graph_algorithm, graph_parameters, varargin)
        %TEST Construct an instance of this class
        %   Detailed explanation goes here

        self.graph_algorithm = graph_algorithm;
        self.graph_parameters = num2cell(graph_parameters);
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'number_of_tests'
                    self.number_of_tests = varargin{2};
                case 'fraction_targets'
                    self.fraction_targets = varargin{2};
                case 'fraction_disturbances'
                    self.fraction_disturbances = varargin{2};
                otherwise
                    error(['Unexpected option: ' varargin{1}])
            end
            varargin(1:2) = [];
        end
    end

    function run(self)
        %METHOD Perform a test according to object properties
        %   Detailed explanation goes here
        
        for i = 1:self.number_of_tests
            G = self.graph_algorithm(self.graph_parameters{:},self.graph_seed + i);
            [self.results_cost(end+1), self.results_time(end+1), self.results_trivial(end+1)] = decouple(G, self.fraction_targets, self.fraction_disturbances);
        end
    end

    function add_to_plot(self)
        % pass
    end
end
end