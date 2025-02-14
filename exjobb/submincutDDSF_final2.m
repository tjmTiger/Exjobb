%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute V_in or V_out (optimal, but not unique) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds V_in or V_out for DDP exploiting the mincut/maxflow algorithm
% Inputs: G connected directed unweighted graph, D disturbance nodes set, T target nodes set
% Additional input: varargin = 'V_in' (finds minimal control nodes), varargin = 'V_out' (finds minimal measured nodes), default = V_in
% Outputs: V_in/V_out optimal in the cardinality sense (not jointly)
% Example : V_in = submincutDDSF_final2(G,D,T,'V_in')
% Additional comment: If I want solutions others than the trivials V_in = T (or V_out = D), I can rewrite T (or D) and G before calling the function
% s.t. T_new (D_new) has K additional in-series nodes with K >= length(V_in) (or length(V_out).
% E.g. A = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 1 1 0], D = 1, T = 4, into A_new = [0 0 0 0 0 0; 1
% 0 0 0 0 0; 1 0 0 0 0 0; 0 1 1 0 0 0; 0 0 0 1 0 0; 0 0 0 1 0 0], T_new = [4 5 6]


function V_in_out = submincutDDSF_final2(G,D,T,varargin)

    if nargin > 3
        a = 2; % for switch
        if strcmp(varargin{1}, 'V_in')
            a = 1;
        elseif strcmp(varargin{1}, 'V_out')
            a = 2;
        else
            disp('Varargin typing error: choose ''V_in'' or ''V_out''');
            return
        end
    else
        a = 1; % V_in as default if varargin empty
    end

    n = numnodes(G);
    k = numedges(G);
    G.Edges.Weight = Inf*ones(k,1);
    n_dist = length(D);
    n_targ = length(T);

    switch a

        case 1
            
            if ~isempty(T) && ~isempty(D)
        
                %%% THIS PART OF THE GRAPH DUPLICANT CONSTRUCTION CAN BE SPEEDED UP %%%

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

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [~,~,cs,ct] = maxflow(G,2*n+1,2*n+2);
                
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
                
                V_in_out = unique(css);
        
            else
                 V_in_out = [];
                 disp('Error: Target or Disturbance set is empty')
                 return
            end

        case 2

            if ~isempty(T) && ~isempty(D)
        
                %%% THIS PART OF THE GRAPH DUPLICANT CONSTRUCTION CAN BE SPEEDED UP %%%

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

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [~,~,cs,ct] = maxflow(G,2*n+1,2*n+2);
                
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
                
                V_in_out = unique(css);
        
            else
                 V_in_out = [];
                 disp('Error: Target or Disturbance set is empty')
                 return
            end

    end    

end