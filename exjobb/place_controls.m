%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute V_in from V_out (suboptimal) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds V_in for DDPOF once V_out is fixed and optimized (by submincutDDSF_final2), the DDP is solved by removing edges only
% Inputs: Adjacency matrix A, Control nodes V_in (computed by submincutDDSF_final2.m), T target nodes set, S_min min conditioned invariant set
% Outputs: V_in input nodes for feedback s.t. im(D) in S_min in Z_max in ker(T)
% Example : V_in = place_controls(A, V_out, T, S_min) 

 function V_in = place_controls(A, V_out, T, S_min) 

     G = digraph(A');
     n = size(A,1);
     p = length(V_out);
     t = length(T);
     V_in = [];
    
     for i = 1:p
         [~,s_in] = findedge(G,outedges(G,V_out(i)));  % find nodes with ingoing edges from V_out
         V_in_tmp = intersect(s_in, setdiff(1:n,S_min)); % Take always an action towards the complement of S_min
         s_in = setdiff(V_in_tmp,V_in); % removes reduntant nodes
         l_tmp = length(s_in);
         for j = 1:l_tmp
             for k = 1:t
                 if ~isempty(shortestpath(G,s_in(j),T(k)))
                     V_in = [V_in; s_in(j)]; % collect such node
                 else
                 end
             end
         end
         V_in = unique(V_in); % remove duplicants
     end
    
     V_in = setdiff(V_in, V_out); % remove output nodes since Vout computed by submincutDDSF_final2.m is s.t. D-T are decoupled, so controlling the out nodes from Vout is enough, no need to cut the edges between the inputs
     V_in = V_in';
     
 end