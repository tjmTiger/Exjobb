%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimizes V_out/V_in conditioned to min(card(V_in))/min(card(V_out)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDOF
% Given the list of optima V_in (or V_out), computes all V_out-V_in
% (V_in-V_out) pairs and outputs the solution that minimizes C = card(V_in)+card(V_out)
% Inputs: G connected directed unweighted graph, D disturbance nodes set, T
% target nodes set, V (V_in or V_out) nodes list of optima solutions computed by mincutDDSF_all.m
% (exhaustive constrained optimization problem if varargin{2} = 'all' in mincutDDSF_all.m)
% Additional inputs: varargin{1} = 'V_in' (given V_out finds V_in), varargin{1} = 'V_out' (given V_in finds V_out)
% Outputs: V_in cell, V_out cell, C cost function (V_in and V_out are cell arrays since the optima may not be unique)
% Example : [V_in, V_out, C] = constrained_optimal_solution(G,D,T,V_out,'V_in'), V_in{i}-V_out{i} is the i-th optimal pair

function [V_in, V_out, C] = constrained_optimal_solution(G,D,T,V,varargin)

    if nargin > 4
        if strcmp(varargin{1}, 'V_in')
            a = 1; % for switch
            % disp('Search for min V_in, given V_out');
            % disp('')
        elseif strcmp(varargin{1}, 'V_out')
            a = 2;
            % disp('Search for min V_out, given V_in');
            % disp('')
        else
            disp('Varargin typing error: choose ''V_in'' or ''V_out'' as additional input');
            return
        end
    else  
        disp('Not enough inputs: choose ''V_in'' or ''V_out'' as additional input');   
        return
    end

    A = full(adjacency(G))';

    switch a

        case 1

            if ~isempty(V)
                
                m = numel(V);
                C = zeros(m,1);
                Cc = C;
                V_in = cell(1, m);

                for i = 1:m
                    Vout = V{i};
                    S_min = ddpf_iff_condition(A, D, T, Vout, 0, 'S_min');
                    Vin = place_controls(A, Vout, T, S_min);
                    V_in{i} = Vin;
                    Cc(i) = length(Vin);
                end

                minValue = min(Cc);
                idx = find(Cc == minValue); % min may be not unique
                C = minValue + length(Vout);
                V_out = V(idx);
                V_in = V_in(idx); % select the optimal subcell

            else
                disp('Error: Target or Disturbance may be empty')
                V_in = [];
                V_out = [];
                C = [];
                return
            end

        case 2

            if ~isempty(V)

                m = numel(V);
                C = zeros(m,1);
                Cc = C;
                V_out = cell(1, m);

                for i = 1:m
                    Vin = V{i};
                    Z_max = ddpf_iff_condition(A, D, T, Vin, 0, 'Z_max');
                    Vout = place_sensors(A, Vin, D, Z_max);
                    V_out{i} = Vout;
                    Cc(i) = length(Vout);
                end

                minValue = min(Cc);
                idx = find(Cc == minValue); % min may be not unique
                C = minValue + length(Vin);
                V_in = V(idx);
                V_out = V_out(idx); % select the optimal subcell

            else
                disp('Error: Target or Disturbance may be empty')
                V_in = [];
                V_out = [];
                C = [];
                return
            end

    end

end