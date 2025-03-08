%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Best constrained optimal V_in-V_out pair when conditioning for V_in or for V_out %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDOF
% Given the list of optimal lists of pairs V_in-V_out provided by
% constrained_optimal_solution.m find min C = card(V_in)+card(V_out), conditioning over V_in or over V_out
% Inputs: [V_in_opt cell, V_out cell, C1 cost function] once constrain
% w.r.t. V_in, and [V_in cell, V_out_opt cell, C2 cost function] once constrain w.r.t. V_out (the order is important)
% Outputs: V_in cell, V_out cell, C cost function, S additional cell telling towards which among V_in or V_out I have conditioned (V_in and V_out are cell arrays since the optima may not be unique)
% S = {'V_in', 'V_in', 'V_out'} meaning that I have 3 optimal V_in-V_out pairs, the first two when conditioning w.r.t. V_in, the last when conditioning w.r.t. V_out
% Example : (the order is important)
% [V_in_opt, V_out, C1] = constrained_optimal_solution(G,D,T,Vout,'V_in');
% [V_in, V_out_opt, C2] = constrained_optimal_solution(G,D,T,Vout,'V_out');
% [V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(V_in_opt, V_out, C1, V_in, V_out_opt, C2)

function [V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(V_in_opt, V_out, C1, V_in, V_out_opt, C2)

    if ~isempty(C1) && ~isempty(C2)
    
        if C1 < C2
            V_in_best = V_in_opt;
            V_out_best = V_out;
            C = C1;
            Ss = ones(numel(V_in_opt),1);
            S = cell(size(Ss)); % initialize the cell array
            S(Ss == 1) = {'V_in'}; % fill 'V_in'
        elseif C2 < C1
            V_in_best = V_in;
            V_out_best = V_out_opt;
            C = C2;
            Ss = ones(numel(V_in),1);
            S = cell(size(Ss)); % initialize the cell array
            S(Ss == 1) = {'V_out'}; % fill 'V_in'
        elseif C1 == C2
            V_in_best = [V_in_opt, V_in];
            V_out_best = [V_out, V_out_opt];
            C = C1;
            x = cell(1,numel(V_in_best));
    
            for i = 1:numel(V_in_best)
                % Sort the arrays in Vin and Vout for comparison
                a = sort(V_in_best{i});
                b = sort(V_out_best{i});
                c = [a b];
                x{i} = c;
            end
    
            V_str = cellfun(@mat2str, x, 'UniformOutput', false);
            [~, idx] = unique(V_str, 'stable');
            allIdx = 1:numel(x);
            repeatedIdx = setdiff(allIdx, idx);
            [~, ~, ic] = unique(V_str, 'stable');
            idxx = find(histcounts(ic, 1:max(ic)+1) > 1);
            repeatedOccurrences = ismember(ic, idxx);
            Ss = [ones(numel(V_in_opt),1); zeros(length(V_in),1)];
            S = cell(size(Ss)); % initialize the cell array
            S(Ss == 1) = {'V_in'}; % fill 'V_in'
            S(Ss == 0) = {'V_out'}; % fill 'V_out'
            S(repeatedOccurrences) = {'V_in and V_out'};
            V_in_best(repeatedIdx) = [];
            V_out_best(repeatedIdx) = [];
            S(repeatedIdx) = [];
        end
    
        T = table(V_in_best', V_out_best', cellfun(@numel, V_in_best)', cellfun(@numel, V_out_best)',  C.*ones(numel(S),1), S, 'VariableNames', {'V_in_best', 'V_out_best','# of Inputs','# of Outputs', 'Minimal Cost C', 'Constraint w.r.t.'});
        disp(T);
    
    else
    
        V_in_best = [];
        V_out_best = [];
        C = [];
        S =[];
        disp('Error: Target or Disturbance may be empty')
    
    end

end