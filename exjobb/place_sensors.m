%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute V_out from V_in (suboptimal) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds V_out for DDPOF once V_in is fixed and optimized (by submincutDDSF_final2), the DDP is solved by removing edges only
% Inputs: Adjacency matrix A, Control nodes V_in (computed by
% submincutDDSF_final2.m), D disturbance nodes set, Z_max max controlled invariant set
% Outputs: V_out output nodes for feedback s.t. im(D) in S_min in Z_max in ker(T)
% Example : V_out = place_sensors(A, V_in, D, Z_max) 

function V_out = place_sensors(A, V_in, D, Z_max)
    
    G = digraph(A');
    m = length(V_in);
    d = length(D);
    V_out = [];
    
    for i = 1:m
        [sout, ~] = findedge(G,inedges(G,V_in(i))); % find nodes with outgoing edges onto V_in
        V_out_tmp = intersect(sout, Z_max); % Take always an action towards Z_max nodes (Necessary? It seems to get (A-GC)V in V)
        sout = setdiff(V_out_tmp,V_out); % removes reduntant nodes
        l_tmp = length(sout);
        for j = 1:l_tmp
            for k = 1:d
                if ~isempty(shortestpath(G,D(k),sout(j)))
                    V_out = [V_out; sout(j)]; % collect such node
                else
                end
            end
        end
        V_out = unique(V_out); % remove duplicants
    end
    
    V_out = setdiff(V_out, V_in); % remove input nodes since Vin computed by submincutDDSF_final2.m is s.t. D-T are decoupled, so observing the in nodes on Vin is enough, no need to cut the edges between the inputs
    V_out = V_out';
    %     [~, Sast, ~] = ddpf_iff_condition(A, D, T, V_in, 0, V_out); % extract only nodes that are in min cond inv superset Sast
    %     V_out = intersect(Sast, V_out);

end