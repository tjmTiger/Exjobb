%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (C,A,B)-pairs construction at minimal cost %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (C,A,B) pairs construction at minimal cost for compensator exploitation for DDP using compensator_ddp.m function
% Inputs: G, D, T
% Outputs: V_out, V_in cell lists for compatible (C,A,B) pairs, s.t.
% exists F|(A-BF)Z_max in Z_max and exists G|(A-GC)S_min in S_min,
% and exists N|(A-BNC)S_min in Z_max, TotalCost (card(V_in)+card(V_out) array), Z_max, S_min
% Example : [V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair(G, D, T);
% Warning: some (C,A,B) pairs could be suitable for DDP after
% precompensation by edge removal through output feedback. This algorithm
% provides only (C,A,B) pairs that do not request precompensation.

function [V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair_backprop_new(G, D, T)
    
    disp('Warning: some (C,A,B) pairs could be suitable for DDP after precompensation by edge removal through output feedback. This algorithm provides only (C,A,B) pairs that do not request precompensation.');
    
    if nargin > 3 || nargin < 3
    
        disp("Error: number of inputs is incorrect, they must be 3 or 4.");
        V_out = [];
        V_in = [];
        TotalCost = [];
        return
    
    end
    
    Vout = submincutDDSF_final2(G,D,T,'V_out');
    A = full(adjacency(G))';
    n = numnodes(G);
    n_dist = length(D);
    n_targ = length(T);
    V_out = [];
    V_in = [];
    TotalCost = [];
    Min_dim = []; % stores intersect(Z,complement(S)) = Z-S = min dim compensator
    S_Min = [];
    Z_Max = [];

    if ~isempty(Vout)
    
        V_out_effective_tot = [];
        G_cell_effective_tot = [];
        S_min_effective_tot = [];
        S_min_out_effective_tot = [];
        V_in_effective_tot = [];
        Z_max_effective_tot = [];
        Vin_cost_tot = [];
    
        %%% V_out_cell_i
        [V_out_cell, G_cell] = mincutDDSF_all(G,D,T,Vout,'V_out','all'); % Provides V_out combinations for DDP at minimum cost
        Vout = []; % save the solutions
        cond1 = 1;
        L2 = 0;

        while cond1 && L2 <= max(length(D),length(T))+5 % second condition in order to faster the search
%         while cond1
    
            xx = zeros(numel(V_out_cell),1);
            L1 = length(V_out_cell{1});
    
            %%% V_out_eff_i
            for i = 1:numel(V_out_cell)
    
                Smin = ddpf_iff_condition(A, D, T, V_out_cell{i}, 0, 'S_min');
                S_min{i} = Smin;
                Smin_out = [];
    
                for ii = 1:length(Smin)
    
                    [~,s_in] = findedge(G,outedges(G,Smin(ii)));
                    Smin_out = [Smin_out; s_in];
    
                end
    
                Smin_out = unique(Smin_out)';
                S_min_out{i} = Smin_out;
    
                %%%%%-----------------------
                if isempty(intersect(Smin_out,T)) % disp('Problem: A*S goes in T and so is not in Zmax')
                    % AAA: This was only necessary if N=0
    
                    xx(i) = 1;
    
                end
                %%%%%-----------------------

                 xx(i) = 1;
    
            end
    
    
            V_out_effective = V_out_cell(xx ~= 0);
            G_cell_effective = G_cell(xx ~= 0);
            S_min_effective = S_min(xx ~= 0);
            S_min_out_effective = S_min_out(xx ~= 0);
%             V_in_effective = cell(numel(V_out_effective),1);
%             Z_max_effective = cell(numel(V_out_effective),1);
%             Vin_cost = zeros(numel(V_out_effective),1);
    
            V_out_effective_tot = [V_out_effective_tot, V_out_effective];
            G_cell_effective_tot = [G_cell_effective_tot, G_cell_effective];
            S_min_effective_tot = [S_min_effective_tot, S_min_effective];
            S_min_out_effective_tot = [S_min_out_effective_tot, S_min_out_effective];
%             V_in_effective_tot = [V_in_effective_tot; V_in_effective];
%             Z_max_effective_tot = [Z_max_effective_tot; Z_max_effective];
%             Vin_cost_tot = [Vin_cost_tot; Vin_cost];
    
            %%% V_out_cell_i+1
            Vout = [Vout, V_out_cell];
            [V_out_cell, G_cell] = mincutDDSF_all(G,D,T,Vout,'V_out','all',G_cell{end});
            L2 = length(V_out_cell{1});
            cond1 = L2-L1 > 0;

            if ~cond1

                result = ones(1, numel(V_out_cell));

                for iii = 1:numel(V_out_cell) % Check if there is any cell in Vout that has the same elements as V_out_cell{i}

                    for jjj = 1:numel(Vout)

                        if length(Vout{jjj}) == length(V_out_cell{iii}) && all(ismember(Vout{jjj}, V_out_cell{iii})) && all(ismember(V_out_cell{iii}, Vout{jjj}))

                            result(iii) = 0;
                            break

                        end

                    end

                end

                if any(result ~= 0)

                    cond1 = 1;
                    nzInd = result ~= 0;
                    V_out_cell = V_out_cell(nzInd);

                end

            end
    
        end
    
        if isempty(V_out_effective_tot)
    
            disp('Problem: the V_out list provided is s.t. A*Smin goes in T and so is not in Zmax. Consider increasing the cost for V_out or add precompensation by edge removal through output feedback.'); % V_out is always contained in S_min
            return
    
        else

            G = digraph(A); % Transpose: invert arrow direction
%             Vin = submincutDDSF_final2(G,T,D,'V_out'); % back-propagation
%             [V_in_cell, G_cell_Vin] = mincutDDSF_all(G,T,D,Vin,'V_out','all'); % back-propagation
    
            for j = 1:numel(V_out_effective_tot) % Allow target and forbid disturbance as control nodes
    
                Smin = S_min_effective_tot{j};
                Smin_out = S_min_out_effective_tot{j};
                
                Vin = submincutDDSF_final2(G,T,D,'V_out'); % back-propagation
                [V_in_cell, G_cell_Vin] = mincutDDSF_all(G,T,D,Vin,'V_out','all'); % back-propagation
    
                cond2 = 1;
                V_in_all = []; % save the solutions

                LL2 = 0;

                while cond2 && LL2 <= max(length(D),length(T))+5 % second condition in order to faster the search

                    %                 while cond2

                    LL1 = length(V_in_cell{1});
                    V_in_all = [V_in_all, V_in_cell];

                    for z = 1:numel(V_in_cell)

                        Zmax = ddpf_iff_condition(A, D, T, V_in_cell{z}, 0, 'Z_max');
                        isSubset1 = all(ismember(Smin, Zmax));
                        isSubset2 = all(ismember(Smin_out, Zmax));
                        isSubset2 = 1;
%                        isSubset3 = isempty(intersect(V_in_cell{z},V_out_effective_tot{j})); % V_in and V_out disjoint
    
%                         if isSubset1 && isSubset2 && isSubset3
                        if isSubset1 && isSubset2
    
                            %%% (C,A,B) pairs
                            V_in_effective = V_in_cell(z);
                            V_out = [V_out; V_out_effective_tot(j)];
                            V_in = [V_in;  V_in_effective];
                            S_Min = [S_Min; {Smin}];
                            Z_Max = [Z_Max; {Zmax}];
                            Total_cost = length(V_out_effective_tot{j})+length(V_in_cell{z});
                            TotalCost = [TotalCost; Total_cost];
                            Min_dim = [Min_dim; length(setdiff(Zmax,Smin))];
    
                        end
    
                    end
    

                    Vin = [Vin, V_in_cell];
                    [V_in_cell, G_cell_Vin] = mincutDDSF_all(G,T,D,Vin,'V_out','all',G_cell_Vin{end});
                    LL2 = length(V_in_cell{1});
                    cond2 = LL2-LL1 > 0;

                    if ~cond2
    
                        result = ones(1, numel(V_in_cell));

                        for iii = 1:numel(V_in_cell) % Check if there is any cell in V_in_all that has the same elements as V_in_cell{i}
    
                            for jjj = 1:numel(V_in_all)
    
                                if length(V_in_all{jjj}) == length(V_in_cell{iii}) && all(ismember(V_in_all{jjj}, V_in_cell{iii})) && all(ismember(V_in_cell{iii}, V_in_all{jjj}))
    
                                    result(iii) = 0;
                                    break
    
                                end
    
                            end
    
                        end
    
                        if any(result ~= 0)
    
                            cond2 = 1;
                            nzInd = result ~= 0;
                            V_in_cell = V_in_cell(nzInd);
    
                        end
    
                    end
    
                end
    
            end
    
        end
    
    else
    
        disp('Problem: V_out = [], DDP already solved, consider to redefine T, D, or G.');
        return
    
    end

    T = table((1:numel(V_out))', V_out, V_in, cellfun(@numel, V_out), cellfun(@numel, V_in),  TotalCost, Min_dim, 'VariableNames', {'(C,A,B)-pair number:', 'V_out_opt', 'V_in_opt', '# of Outputs', '# of Inputs', 'Minimal Cost C','Minimal reduced compensator dimension'});
    disp(T);
    
end