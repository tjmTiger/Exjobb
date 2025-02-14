%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTROL INVARIANCE IFF CONDITION FOR DDPS/OF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IFF CONDITION FOR DDPSF: Construct Zast = max controlled invariant set in ker(T)
% IFF CONDITION FOR DDPOF: If varargin = V_out: Construct Sast = min conditioned invariant set containing V_dist
% Can also be exploited to compute Zmax or Sast only
% Inputs: Adjacency matrix A, Distrurbance nodes V_dist, Target nodes V_targ, Input nodes V_in, show variable (if 1 display, otherwise not)
% Additional input: varargin{1} = Output nodes V_out (DDPOF), varargin{1} = 'Z_max' compute only Z_max, varargin{1} = 'S_min' compute only S_min
% Outputs: Zast max control invariant subset, Sast min conditioned invariant superset (if varargin not empty, otherwise Sast = []),
% d = Boolean variable if V_out is empty otherwise d = [d t] with d,t Boolean variables
% (d = 1 if DDPSF solvable, 0 otherwise, t = 1 if DDPOF solvable, 0 otherwise).
% The function can give the three outputs [Zast, Sast, d] or only one of them Zast or Sast if requested.
% Example 1: [Zast, ~, d] = ddpf_iff_condition(A, V_dist, V_targ, V_in, 0) % no show, DDPSF (Sast = [])
% Example 2: Sast = ddpf_iff_condition(A, V_dist, V_targ, V, 0, 'S_min') % no show, compute S_min, V = V_out

function varargout = ddpf_iff_condition(A, V_dist, V_targ, V, show, varargin)
        
    if nargin > 5
        if strcmp(varargin{1}, 'Z_max')
            V_in = V;
            a = 3;
        elseif strcmp(varargin{1}, 'S_min')
            V_out = V;
            a = 2;
        else
            V_in = V;
            V_out = varargin{1};
            a = 1;
        end
    else
        V_in = V;
        V_out = [];
        a = 1;
    end
    
    switch a
    
        case 1
    
            if isempty(V_out)
                Sast = [];
                disp('V_out is empty')
            end
    
            n = size(A,1); % number of nodes
            Zo = setdiff(1:n,V_targ); % remove from total nodes the target nodes
            Zcurr = Zo;
            flag_done = 0;
    
            while ~flag_done
                Zelim = [];
                for i = 1:length(Zcurr)
                    ind = find(A(:,Zcurr(i)));
                    if ~isempty(ind)
                        flag = ismember(ind,union(Zcurr,V_in)); % check if all heads are in Zcurr or in V_in
                        if sum(flag) < length(ind) % some are not
                            Zelim = union(Zelim, Zcurr(i));
                        end
                    end
                end
                if isempty(Zelim)
                    flag_done = 1;
                    Zast = Zcurr;
                    if show == 1
                        disp(['Size of Zast = ',num2str(length(Zcurr))]);
                    end
                else
                    if show == 1
                        disp(['Elimininate ',num2str(length(Zelim)),' out of ',num2str(length(Zcurr))]);
                    end
                    Zcurr = setdiff(Zcurr,Zelim);
                end
            end
    
            % check DDP by state feedback
    
            flag_D_in_Zast = ismember(V_dist,Zast);
    
            if sum(flag_D_in_Zast) == length(V_dist)
                d = 1;
                if show == 1
                    disp('The DDPSF problem is solvable');
                end
            else
                d = 0;
                if show == 1
                    disp(['The DDPSF problem is not solvable: ', num2str(sum(flag_D_in_Zast)),' out of ',...
                        num2str(length(V_dist)),' are decouplable with these ',num2str(length(V_in)),' controls']);
                end
            end
    
            if ~isempty(V_out)
                So = V_dist;
                Scurr = So;
                flag_done = 0;
                while ~flag_done
                    Sadd = [];
                    Scurr_no_out = setdiff(Scurr,V_out);
                    if ~isempty(Scurr_no_out)
                        for i = 1:length(Scurr_no_out)
                            ind = find(A(:,Scurr_no_out(i))); % find all heads of edges with tail in Scurr-V_out
                            flag = ismember(ind,Scurr); % check if all heads are in Scurr
                            if sum(flag) < length(ind) % some are not
                                Sadd = union(Sadd, setdiff(ind,Scurr));
                            end
                        end
                    end
                    if isempty(Sadd)
                        flag_done = 1;
                        Sast = Scurr;
                        if show == 1
                            disp(['Size of Sast = ',num2str(length(Scurr))])
                        end
                    else
                        if show == 1
                            disp(['Add ',num2str(length(Sadd)),' to Scurr of size ',num2str(length(Scurr))])
                        end
                        Scurr = union(Scurr,Sadd);
                    end
                end
    
                % check DDP by output feedback
    
                flag_Sast_in_Zast = ismember(Sast,Zast);
                if sum(flag_Sast_in_Zast) == length(Sast)
                    t = 1;
                    if show == 1
                        disp('The DDPOF problem is solvable by output feedback')
                    end
                else
                    t = 0;
                    if show == 1
                        disp('The DDPOF problem is not solvable by output feedback: ')
                        disp(['dim(Sast)= ',num2str(length(Sast)),', dim(Zast)= ',num2str(sum(flag_Sast_in_Zast)),...
                            ', with ',num2str(sum(flag_Sast_in_Zast)),' nodes of Sast in Zast, ',...
                            'with ',num2str(length(V_in)),' controls'])
                    end
                end
                d = [d t];
            end
    
        case 2
    
            if ~isempty(V_out)
                So = V_dist;
                Scurr = So;
                flag_done = 0;
                while ~flag_done
                    Sadd = [];
                    Scurr_no_out = setdiff(Scurr,V_out);
                    if ~isempty(Scurr_no_out)
                        for i = 1:length(Scurr_no_out)
                            ind = find(A(:,Scurr_no_out(i))); % find all heads of edges with tail in Scurr-V_out
                            flag = ismember(ind,Scurr); % check if all heads are in Scurr
                            if sum(flag) < length(ind) % some are not
                                Sadd = union(Sadd, setdiff(ind,Scurr));
                            end
                        end
                    end
                    if isempty(Sadd)
                        flag_done = 1;
                        Sast = Scurr;
                        if show == 1
                            disp(['Size of Sast = ',num2str(length(Scurr))])
                        end
                    else
                        if show == 1
                            disp(['Add ',num2str(length(Sadd)),' to Scurr of size ',num2str(length(Scurr))])
                        end
                        Scurr = union(Scurr,Sadd);
                    end
                end
            else
                disp('Error: V_out is empty')
            end
    
        case 3
    
            if ~isempty(V_in)
                n = size(A,1); % number of nodes
                Zo = setdiff(1:n,V_targ); % remove from total nodes the target nodes
                Zcurr = Zo;
                flag_done = 0;
    
                while ~flag_done
                    Zelim = [];
                    for i = 1:length(Zcurr)
                        ind = find(A(:,Zcurr(i)));
                        if ~isempty(ind)
                            flag = ismember(ind,union(Zcurr,V_in)); % check if all heads are in Zcurr or in V_in
                            if sum(flag) < length(ind) % some are not
                                Zelim = union(Zelim, Zcurr(i));
                            end
                        end
                    end
                    if isempty(Zelim)
                        flag_done = 1;
                        Zast = Zcurr;
                        if show == 1
                            disp(['Size of Zast = ',num2str(length(Zcurr))]);
                        end
                    else
                        if show == 1
                            disp(['Elimininate ',num2str(length(Zelim)),' out of ',num2str(length(Zcurr))]);
                        end
                        Zcurr = setdiff(Zcurr,Zelim);
                    end
                end
    
            else
                disp('Error: V_in is empty')
            end
    end
    
    nargoutchk(1, 3);  % check that the number of outputs is between 1 and 3
    
    switch nargout % assign outputs based on the number requested

        case 3

            varargout{1} = Zast;
            varargout{2} = Sast';
            varargout{3} = d;
%             return;

        case 1

            if nargin > 5
                if strcmp(varargin{1}, 'Z_max')
                    varargout{1} = Zast;
                elseif strcmp(varargin{1}, 'S_min')
                    varargout{1} = Sast';
                else
                    disp('Varargin typing error: choose ''Z_max'' or ''S_min'' as first additional input')
                    return
                end
            else
                disp('Error: additional input missing, assign ''Z_max'' or ''S_min'' in varargin{1}')
                return
            end

    end

end