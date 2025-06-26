classdef Tests < handle
%TESTS Summary of this class goes here
%   Detailed explanation goes here

properties
    % number_of_tests {mustBeNumeric} = 200 => fix so that all tasts have
    % same sample size !!!!
    tests {mustBeCell} = {};
end

methods
    function self = Tests()
        %TESTS Construct an instance of this class
        %   Detailed explanation goes here
    end
    
    function run(self,varargin)
        %METHOD1 Summary of this method goes here
        %   Detailed explanation goes here
        t = Test(varargin{:});
        t.run()
        self.tests{end+1} = t;
    end

    function plot(self)
        %METHOD2 Summary of this method goes here
        %   Detailed explanation goes here
        disp("todo: implement plotting")
    end
end
end

