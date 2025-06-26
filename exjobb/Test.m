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
    sample_size {mustBeNumeric} = 200
    fraction_targets {mustBeFloat} = 0.1
    fraction_disturbances {mustBeFloat} = 0.1

    % other
    results_cost
    results_time
    results_trivial
end
methods
    function self = Test(graph_algorithm, graph_parameters, varargin)
        %TEST Construct an instance of this class
        %   Detailed explanation goes here

        self.graph_algorithm = graph_algorithm;
        self.graph_parameters = num2cell(graph_parameters);
        while ~isempty(varargin)
            switch lower(varargin{1})
                case 'sample_size'
                    self.sample_size = varargin{2};
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

        % parfor does not like self, change all variables we use to local
        graph_algorithm = self.graph_algorithm;
        graph_parameters = self.graph_parameters;
        graph_seed = self.graph_seed;
        results_cost = zeros(1,self.sample_size);
        results_time = zeros(1,self.sample_size);
        results_trivial = zeros(1,self.sample_size);
        parfor i = 1:self.sample_size
            G = graph_algorithm(graph_parameters{:},graph_seed + i);
            [results_cost(i), results_time(i), results_trivial(i)] = decouple(G, self.fraction_targets, self.fraction_disturbances);
        end
        self.results_cost = results_cost;
        self.results_time = results_time;
        self.results_trivial = results_trivial;
    end

    function add_to_plot(self)
        % pass
    end
end
end