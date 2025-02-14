%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTROL INVARIANCE FRIEND COMPUTATION FOR DDPS/OF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FRIEND F COMPUTATION: Makes the new state matrix A-BF invariant (given matrix V with column as basis of Z_max (max control invariant subset,
% the matrix inv(T)*A*T is of the form [A11 A12; 0 A22] with T=[V X] and X complement of the basis V)
% Inputs: Adjacency matrix A, Control nodes V_in (computed by submincutDDSF_final2.m),
% Z_max max control invariant subset (computed by ddpf_iff_condition.m), D disturbance nodes set, T target nodes set
% Additional input: varargin = Output nodes V_out (DDPOF) (computed by place_sensors.m) s.t. An-G*Cn invariant by output injection with friend G
% Outputs: F, An, Bn, Tn, Dn new state-space matrices after change of basis, TT chane of basis matrix
% Example : [F, An, Bn, Tn, Dn, TT, Cn, G, Cn] = compute_friend(A, V_in, Z_max, D, T, V_out) 

function [F, An, Bn, Tn, Dn, TT, G, Cn] = compute_friend(A, V_in, Z_max, D, T, varargin)

    if nargin > 5
        a = 2; % for switch, output injection
        V_out = varargin{1};
        y = length(V_out);
    else
        a = 1; % state feedback
        V_out = [];
        G = [];
        Cn = [];
    end
    
    switch a
    
        case 1
    
            % Compute change of basis to have a left-invariant system (s.t. pseudoinverse = transpose)
            n = size(A,1);
            s = length(Z_max);
            m = length(V_in);
            d = length(D);
            t = length(T);
            K = unique([Z_max V_in],'stable'); % [Z_max V_in];
%           c = s+m;
            c = length(K);
            TT = zeros(n,n);
    
            for i = 1:c
                TT(K(i),i) = 1;
            end
    
            P = setdiff(1:1:n,K); % Complementary set
    
            for i = 1:(n-c)
                TT(P(i),i+c) = 1;
            end
    
            Bi = zeros(n,m); % Initial B matrix
    
            for i = 1:m
                Bi(V_in(i),i) = 1;
            end
    
            Di = zeros(n,d); % Initial D matrix
    
            for i = 1:d
                Di(D(i),i) = 1;
            end
    
            Ti = zeros(t,n); % Initial T matrix
    
            for i = 1:t
                Ti(i,T(i)) = 1;
            end
    
            An = inv(TT)*A*TT;
            Bn = inv(TT)*Bi;
            Dn = inv(TT)*Di;
            Tn = Ti*TT;
    
            %     V = zeros(n,s); % Compute basis for Z_max
            %
            %     for i = 1:s
            %         V(Z_max(i),i) = 1;
            %     end
    
            V = [eye(s,s);zeros(n-s,s)]; % Basis for Z_max with new relabeling
            Kn = [V Bn];
    
            AV = An*V;
            X = Kn'*AV; % X = pinv(K)*AV;
            U = X(end-(m-1):end,:);
            F = U*V'; % F = U*inv(V'*V)*V'; % Friend is unique since system left-invariant
    
        case 2
    
            % Compute change of basis to have a left-invariant system (s.t. pseudoinverse = transpose)
            n = size(A,1);
            s = length(Z_max);
            m = length(V_in);
            d = length(D);
            t = length(T);
            K = unique([Z_max V_in],'stable'); % [Z_max V_in];
%           c = s+m;
            c = length(K);
            TT = zeros(n,n);
    
            for i = 1:c
                TT(K(i),i) = 1;
            end
    
            P = setdiff(1:1:n,K); % Complementary set
    
            for i = 1:(n-c)
                TT(P(i),i+c) = 1;
            end
    
            Bi = zeros(n,m); % Initial B matrix
    
            for i = 1:m
                Bi(V_in(i),i) = 1;
            end
    
            Di = zeros(n,d); % Initial D matrix
    
            for i = 1:d
                Di(D(i),i) = 1;
            end
    
            Ti = zeros(t,n); % Initial T matrix
    
            for i = 1:t
                Ti(i,T(i)) = 1;
            end
    
            Ci = zeros(y,n); % Initial C matrix
    
            for i = 1:y
                Ci(i,V_out(i)) = 1;
            end
    
            An = inv(TT)*A*TT;
            Bn = inv(TT)*Bi;
            Dn = inv(TT)*Di;
            Tn = Ti*TT;
            Cn = Ci*TT;
    
            V = [eye(s,s); zeros(n-s,s)]; % Basis for Z_max with new relabeling
            Kn = [V Bn];
    
            AV = An*V;
            X = Kn'*AV; % X = pinv(Kn)*AV;
            U = X(end-(m-1):end,:);
            F = U*V'; % F = U*inv(V'*V)*V'; % Friend is unique since system left-invariant
    
            [~, col] = find(Cn); % new positions of V_out nodes after change of basis
            G = zeros(n,y);
            [V_in_new, ~] = find(Bn);
            G_new = digraph(An');
    
            for i = 1:m
                [sout, ~] = findedge(G_new,inedges(G_new,V_in_new(i))); % find nodes with outgoing edges onto V_in
                sout = setdiff(sout,V_in_new); % neglect self-loops
                col_i = intersect(col,sout);
                for j = 1:y
                    [~, col_j] = find(Cn(j,:));
                    if ~isempty(intersect(col_j,col_i))
                        G(V_in_new(i),j) = An(V_in_new(i),col_j);
                    else
                        G(V_in_new(i),j) = 0;
                    end
                end
            end
    
    end

end