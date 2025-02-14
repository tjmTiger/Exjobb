%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute all optimal V_in or V_out (optimal, but not unique) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds all optima V_in or V_out for DDP exploiting the mincut/maxflow algorithm
% Inputs: G connected directed unweighted graph, D disturbance nodes set, T
% target nodes set, V_in/V_out nodes optimal solution computed by submincutDDSF_final2.m
% Additional inputs: varargin{1} = 'V_in' (finds all minimal control nodes), varargin{1} = 'V_out' (finds all minimal measured nodes), default = V_in
% varargin{2} = 'all' or 'partial' (finds the set of solutions, 'all' is a
% combinatorial NP-hard search, while 'partial' has a computational cost of M*I*N,
% where M is the length of the solution, I is the maxflow-mincut
% polynomial time solution cost, and N is the number of optimal solutions)
% varargin{3} = G_extended the algorithm searches all solutions equi-cost
% by starting with the preset weighted graph G_extended (ex. G_extended =
% G_cell{end}, and V_in_out = V{end} outputs of a previous run of mincutDDSF_all.m)
% Outputs: V_in/V_out optimal in the cardinality sense (not jointly) and G_cell cell array of all weighted graphs for all optimal solutions
% Example : V_in_all = mincutDDSF_all(G,D,T,V_in,'V_in','partial')
% Additional comment: If I want solutions others than the trivials V_in = T (or V_out = D), I can rewrite T (or D) and G before calling the function
% s.t. T_new (D_new) has K additional in-series nodes with K >= length(V_in) (or length(V_out).
% E.g. A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0], D = 1, T = 4, into A_new = [0 0 0 0 0 0; 1
% 0 0 0 0 0; 1 0 0 0 0 0; 0 1 1 0 0 0; 0 0 0 1 0 0; 0 0 0 1 0 0], T_new = [4 5 6]
% Warning: when the function has 7 inputs (i.e. varargin{3} = G_extended,
% then also the 4th input V_in_out must be a cell array compatible with the cell G_extended)


function varargout = mincutDDSF_all(G,D,T,V_in_out,varargin)

    if nargin > 4
        if strcmp(varargin{1}, 'V_in')
            a = 1; % for switch
            disp('Search for min V_in');
            disp('')
        elseif strcmp(varargin{1}, 'V_out')
            a = 2;
            disp('Search for min V_out');
            disp('')
        else
            disp('Varargin typing error: choose ''V_in'' or ''V_out'' as first additional input and ''all'' or ''partial'' as second additional input');
            return
        end
        if nargin > 5
            if strcmp(varargin{2}, 'all')
                b = 2; % for switch
                disp('Search for ''all'' optimal solutions');
                disp('')
            elseif strcmp(varargin{2}, 'partial')
                b = 1;
                disp('Search for some optimal solutions');
                disp('')
            else
                disp('Varargin typing error: choose ''V_in'' or ''V_out'' as first additional input and ''all'' or ''partial'' as second additional input');
                return
            end
        else
            b = 1;
        end
        if nargin == 7
            a = 3; % use preset G_extended structure
            G_extended = varargin{3};
            b = 3; % search for all suboptimal solutions starting by G_extended, maxflow value is increased
            disp('Warning: the solutions are found by starting with a preset weighted G_extended graph. Such solutions may not be optimal in the sense of the cardinality.')
        elseif nargin > 7
            disp('Error: too many inputs. Max nargin = 7.')
            return
        end
    elseif nargin == 4
        a = 1; % 'V_in' as default if varargin empty
        b = 1; % 'partial' as default if varargin empty
    else
        disp('Not enough inputs: choose ''V_in'' or ''V_out'' as first additional input and ''all'' or ''partial'' as second additional input');
        return
    end

    n = numnodes(G);
    k = numedges(G);
    G.Edges.Weight = Inf*ones(k,1);
    n_dist = length(D);
    n_targ = length(T);
    
    switch a
    
        case 1
    
            if ~isempty(T) && ~isempty(D)
    
                for i = 1:n
                    % G = addnode(G,i+n);
                    [~,nid] = outedges(G,i);
                    for j = 1:length(nid)
                        G = rmedge(G,i,nid(j));
                        G = addedge(G,i+n,nid(j),Inf);
                    end
                    if isempty(intersect(i,D))
                        G = addedge(G,i,i+n,1);
                    else
                        G = addedge(G,i,i+n,Inf);
                    end
                end
    
                for i = 1:n_dist
                    G = addedge(G,2*n+1,D(i),Inf);
                end
    
                for i = 1:n_targ
                    G = addedge(G,T(i)+n,2*n+2,Inf);
                end
    
            else
                V = {};
                disp('Error: Target or Disturbance set is empty')
                return
            end
    
        case 2
    
            if ~isempty(T) && ~isempty(D)
    
                for i = 1:n
                    % G = addnode(G,i+n);
                    [~,nid] = outedges(G,i);
                    for j = 1:length(nid)
                        G = rmedge(G,i,nid(j));
                        G = addedge(G,i+n,nid(j),Inf);
                    end
                    if isempty(intersect(i,T))
                        G = addedge(G,i,i+n,1);
                    else
                        G = addedge(G,i,i+n,Inf);
                    end
                end
    
                for i = 1:n_dist
                    G = addedge(G,2*n+1,D(i),Inf);
                end
    
                for i = 1:n_targ
                    G = addedge(G,T(i)+n,2*n+2,Inf);
                end
    
            else
                V = {};
                disp('Error: Target or Disturbance set is empty')
                return
            end
    
        case 3
    
            G = G_extended;
    
    end
    
    switch b
    
        case 1
    
            maxFlow = length(V_in_out); % works if V_in_out is optimal solution computed via submincutDDSF_final2.m
            maxFlow_n = maxFlow;
            ii = 1;
            V{ii} = V_in_out; % collect solutions in a cell
            mu = 0.5; % parameter < 1 strictly
            p = mu/maxFlow;
    
            while maxFlow_n == maxFlow
    
                for i = 1:maxFlow
                    G.Edges.Weight(V_in_out(i)) = 1+p;
                end
    
                [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
    
                if maxFlow_n > maxFlow
                    if a == 1
                        disp('V_in is the unique optimal solution provided by the algorithm')
                        disp('')
                    elseif a == 2
                        disp('V_out is the unique optimal solution provided by the algorithm')
                        disp('')
                    end
                    return
                end
    
                if ~ismember(2*n+1,cs) % make sure that cs set is from the disturbance side
                    ck = cs;
                    cs = ct;
                    ct = ck;
                end
    
                css = [];
                ctt = [];
    
                for i = 1:length(cs)
                    for j = 1:length(ct)
                        if findedge(G,cs(i),ct(j)) > 0
                            css = [css, cs(i)];
                            ctt = [ctt, ct(j)];
                        end
                    end
                end
    
                V_in_out_new = unique(css);
                V_in_out = V_in_out_new;
                ii = ii+1;
                V{ii} = V_in_out;
    
            end
    
        case 2
    
            G_cell{1} = G;
            maxFlow = length(V_in_out); % works if V_in_out is optimal solution computed via submincutDDSF_final2.m
            maxFlow_n = maxFlow;
            ii = 1;
            kk = ii;
            V_is_increasing = 1;
            V{ii} = V_in_out; % collect solutions in a cell
            mu = 0.5; % parameter < 1 strictly
            p = mu/maxFlow;
            MMaxFlow = [];
    
            while V_is_increasing
    
                V_in_out = V{kk};
    
                for s = 1:maxFlow
    
                    G.Edges.Weight(V_in_out(s)) = 1+p; % inf;
                    [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                    MMaxFlow = [MMaxFlow; maxFlow_n];
    
                    if maxFlow_n == maxFlow
                        if ~ismember(2*n+1,cs) % make sure that cs set is from the disturbance side
                            ck = cs;
                            cs = ct;
                            ct = ck;
                        end
    
                        css = [];
                        ctt = [];
    
                        for i = 1:length(cs)
                            for j = 1:length(ct)
                                if findedge(G,cs(i),ct(j)) > 0
                                    css = [css, cs(i)];
                                    ctt = [ctt, ct(j)];
                                end
                            end
                        end
    
                        V_in_out_new = unique(css);
                        ii = ii+1;
                        V{ii} = V_in_out_new;
                        G_cell{ii} = G;
    
                    end
    
                    G.Edges.Weight(V_in_out(s)) = 1; % back in place
    
                end
    
                %                 if all(MMaxFlow ~= maxFlow)
                %                     V_is_increasing = 0;
                %                     return % finish
                %                 end
    
                for i = 1:maxFlow
                    G.Edges.Weight(V_in_out(i)) = 1+p;
                end
    
                [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                MMaxFlow = [MMaxFlow; maxFlow_n];
    
                if maxFlow_n == maxFlow
                    if ~ismember(2*n+1,cs) % make sure that cs set is from the disturbance side
                        ck = cs;
                        cs = ct;
                        ct = ck;
                    end
    
                    css = [];
                    ctt = [];
    
                    for i = 1:length(cs)
                        for j = 1:length(ct)
                            if findedge(G,cs(i),ct(j)) > 0
                                css = [css, cs(i)];
                                ctt = [ctt, ct(j)];
                            end
                        end
                    end
    
                    V_in_out_new = unique(css);
                    ii = ii+1;
                    V{ii} = V_in_out_new;
                    G_cell{ii} = G;
    
                end
    
                V_sor = cellfun(@sort, V, 'UniformOutput', false);
                V_str = cellfun(@mat2str, V_sor, 'UniformOutput', false);
                [~, idx] = unique(V_str, 'stable');
                V = V(idx);
                V = V(~cellfun('isempty', V));
                G_cell = G_cell(idx);
                G_cell = G_cell(~cellfun('isempty', G_cell));
                kk = kk+1;
    
                if kk > numel(V) % no more solutions with maxFlow-cost are found
                    varargout{1} = V;
                    varargout{2} = G_cell;
                    return
                end
    
            end
    
        case 3

            maxFlow_in = length(V_in_out{end});
            mu_preset = 1.1;
            VV = [];

            for i = 1:numel(V_in_out)
                V_tmp = V_in_out{i};
                L = length(V_tmp);
                p_preset = mu_preset/L;
                for j = 1:L
                     if isempty(intersect(V_tmp(j),VV)) % updated weight added only once to each node in previous solutions sets
                        G.Edges.Weight(V_tmp(j)) = G.Edges.Weight(V_tmp(j))+p_preset;
                     end
                end
                VV = [V_tmp, VV];
            end

            [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
            mu_increase = 1; 
            condition = 1;
            tol = 0.09;
            k = 1;
            par = 2;

            while condition
                VV = []; % reset
                if abs(1+maxFlow_in-maxFlow_n) <= tol
                    condition = 0;
                    maxFlow = maxFlow_n;
                else
                    for i = 1:numel(V_in_out)
                        V_tmp = V_in_out{i};
                        L = length(V_tmp);
                        p_increase = mu_increase/L;
                        for j = 1:L
                            if isempty(intersect(V_tmp(j),VV)) % updated weight added only once to each node in previous solutions sets
                                G.Edges.Weight(V_tmp(j)) = G.Edges.Weight(V_tmp(j))+p_increase;
                            end
                        end
                        VV = [V_tmp, VV];
                    end
                    [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                    maxFlow_in = maxFlow_in + 1;
                end
                k = k + 1;
                if k == par*n
                    disp('Warning: unable to compute new solutions. Probably the solution already exists.')
                    maxFlow = maxFlow_n;
%                     varargout{1} = [];
%                     varargout{2} = [];
                    break
                end
            end

            if ~ismember(2*n+1,cs) % make sure that cs set is from the disturbance side
                ck = cs;
                cs = ct;
                ct = ck;
            end

            css = [];
            ctt = [];

            for i = 1:length(cs)
                for j = 1:length(ct)
                    if findedge(G,cs(i),ct(j)) > 0
                        css = [css, cs(i)];
                        ctt = [ctt, ct(j)];
                    end
                end
            end

            V_in_out_new = unique(css);
            V{1} = V_in_out_new; % substitute initial V_in_out with new solutions of cost maxFlow > maxFlow_in since ii initialized at 0
            G_cell{1} = G;
            ii = 1;
            kk = 1;
            V_is_increasing = 1;
            mu = 0.5; % parameter < 1 strictly
            p = mu/maxFlow;
            MMaxFlow = [];
    
            while V_is_increasing
    
                V_in_out = V{kk};
    
                for s = 1:length(V_in_out)
    
                    G.Edges.Weight(V_in_out(s)) = G.Edges.Weight(V_in_out(s))+p; % inf;
                    [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                    MMaxFlow = [MMaxFlow; maxFlow_n];
    
                    if maxFlow_n == maxFlow
                        if ~ismember(2*n+1,cs) % make sure that cs set is from the disturbance side
                            ck = cs;
                            cs = ct;
                            ct = ck;
                        end
    
                        css = [];
                        ctt = [];
    
                        for i = 1:length(cs)
                            for j = 1:length(ct)
                                if findedge(G,cs(i),ct(j)) > 0
                                    css = [css, cs(i)];
                                    ctt = [ctt, ct(j)];
                                end
                            end
                        end
    
                        V_in_out_new = unique(css);
                        ii = ii+1;
                        V{ii} = V_in_out_new; 
                        G_cell{ii} = G;
    
                    end
    
                    G.Edges.Weight(V_in_out(s)) = G.Edges.Weight(V_in_out(s))-p; % back in place
    
                end
    
                %                 if all(MMaxFlow ~= maxFlow)
                %                     V_is_increasing = 0;
                %                     return % finish
                %                 end
    
                for i = 1:length(V_in_out)
                    G.Edges.Weight(V_in_out(i)) = G.Edges.Weight(V_in_out(i))+p;
                end
    
                [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                MMaxFlow = [MMaxFlow; maxFlow_n];
    
                if maxFlow_n == maxFlow
                    if ~ismember(2*n+1,cs) % make sure that cs set is from the disturbance side
                        ck = cs;
                        cs = ct;
                        ct = ck;
                    end
    
                    css = [];
                    ctt = [];
    
                    for i = 1:length(cs)
                        for j = 1:length(ct)
                            if findedge(G,cs(i),ct(j)) > 0
                                css = [css, cs(i)];
                                ctt = [ctt, ct(j)];
                            end
                        end
                    end
    
                    V_in_out_new = unique(css);
                    ii = ii+1;
                    V{ii} = V_in_out_new;
                    G_cell{ii} = G;
    
                end
    
                V_sor = cellfun(@sort, V, 'UniformOutput', false);
                V_str = cellfun(@mat2str, V_sor, 'UniformOutput', false);
                [~, idx] = unique(V_str, 'stable');
                V = V(idx);
                V = V(~cellfun('isempty', V));
                G_cell = G_cell(idx);
                G_cell = G_cell(~cellfun('isempty', G_cell));
                kk = kk+1;
    
                if kk > numel(V) % no more solutions with maxFlow-cost are found
                    varargout{1} = V;
                    varargout{2} = G_cell;
                    return
                end
    
            end
    
    end   
    
    nargoutchk(1, 2);  % check that the number of outputs is between 1 and 3
    
    switch nargout % assign outputs based on the number requested
    
        case 2
    
            if b == 2
                varargout{1} = V;
                varargout{2} = G_cell;
            else
                disp('The G_cell is displayed only if varargin{2} = ''all''')
                return;
            end
    
        case 1
    
            varargout{1} = V;
    
    end

end
