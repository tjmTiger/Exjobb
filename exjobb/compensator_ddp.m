%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compensator synthesis using (C,A,B)-pairs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (full state)-COMPENSATOR synthesis with minimal input-ouput (C,A,B) pair: s.t. d/dt[x;w] = A_e*[x;w]+D_e*d,
% t = T_e*[x;w] and T_e*exp(A_e*t)*D_e = 0, forall t, or equiv.T_e*inv(sI-A_e)*D_e = 0, forall s. A_e = [A+BNC BM; LC K], D_e = [D;0], T_e = [T 0].
% If exists F|(A-BF)Z_max in Z_max and exists G|(A-GC)S_min in S_min and exists N|(A-BNC)S_min in Z_max, then full state observer with
% K = A-BF-GC+BNC, L = G-BN, M = NC -F. The proposed algorithm finds optimal B,C s.t. S_min in Z_max.
% Inputs: G, D, T, V_out, V_in, S_Min, Z_Max (last 4 inputs from cab_pair.m)
% Outputs: COM compensator structure COM.A_e, COM.D_e, COM.T_e, COM.V_in, COM.V_out, COM.GAMMA (gamma = [K L; M N] compensator) for different (C,A,B) pairs
% Example : COM = compensator_ddp(G, D, T, V_out, V_in, S_Min, Z_Max)
% Additional comment: The algorithm needs compatible (C,A,B)-pairs (consider precompute using cab_pair.m)

function COM = compensator_ddp(G, D, T, V_out, V_in, S_Min, Z_Max)

    if ~isempty(V_out)
    
        A = full(adjacency(G))';
        Vout = place_sensors(A, V_in, D, Z_Max);
        [F, An, Bn, Tn, Dn, TT1, ~, ~] = compute_friend(A, V_in, Z_Max, D, T);
        % Transform back in the original basis
        Ao1 = TT1*An*inv(TT1); % should be equal to A
        Bo = TT1*Bn;
        To = Tn*inv(TT1);
        Do = TT1*Dn;
        Fo = F*inv(TT1);
        Vin = place_controls(A, V_out, T, S_Min);
        Zmax = ddpf_iff_condition(A, D, T, Vin, 0, 'Z_max');
        [~, ~, ~, ~, ~, TT2, Gg, Cn] = compute_friend(A, Vin, Zmax, D, T, V_out);
        %         [~, An2, ~, ~, ~, TT2, Gg, Cn] = compute_friend(A, Vin, Zmax, D, T, V_out);
    
        N = zeros(length(V_in), length(V_out));
        for i = 1:length(V_in)
            for j = 1:length(V_out)
                if ~isempty(A(V_in(i),V_out(j)))
                    N(i,j) = A(V_in(i),V_out(j));
                end
            end
        end
        % Transform back in the original basis
        %         Ao2 = TT2*An2*inv(TT2); % should be equal to A
        Co = Cn*inv(TT2);
        Go = TT2*Gg;
    
        %         AA1 = Ao1-Bo*Fo;
        %         AA2 = Ao2-Go*Co;
        %
%         figure
%         plot(digraph(Ao1'))
%     
%         figure
%         plot(digraph(Ao2'))
%     
%         figure
%         plot(digraph(AA1'))
%     
%         figure
%         plot(digraph(AA2'))
%     
%         Ao = A;
    
%         N = zeros(length(V_in), length(V_out));
        K = A-Bo*Fo-Go*Co+Bo*N*Co;
        L = -Bo*N+Go;
        M = -Fo+N*Co;
    
        gamma = [K L; M N];
        Ae = [A-Bo*N*Co Bo*M; L*Co K];
        De = [Do; zeros(size(Do,1),size(Do,2))];
        Te = [To zeros(size(To,1),size(To,2))];
    
        COM.A_e = Ae; COM.D_e = De; COM.T_e = Te; COM.V_in = V_in; COM.V_out = V_out; COM.GAMMA = gamma;
    
    else
        disp('Problem: V_out = [], DDP already solved, consider to redefine T, D, or G.')
        COM = [];
        return
    end

end