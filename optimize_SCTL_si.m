% optimize_SCTL_si estimates the random effects parameters and computes
% derivatives wrt beta and delta moreover some derivatives of the hessian
% of the respective objective function are also returned as they are needed
% for certain computations later on and computational complexity is rather
% high
%
% USAGE:
% ======
% [...] = bhat_SCTL_si(Model,Data,bhat_0,beta,delta,type_D,t,Ym,Tm,ind_y,ind_t,F_diff,b_diff,s,fms)
%
% INPUTS:
% =======
% Model ... [struct] model definition
% Data ... [struct] data definition
% bhat_0 ... [1xnb] previously estimated random effect parameters, these are used
%     as initialistion to the next estimation
% beta ... [1xnbeta] common effect parameter
% delta ... [1xndelta] parametrisation of random effect covariance
% type_D ... [string] covariance parametrisation type
% t ... [1xnt]time vector for simulation
% Ym ... [ntxnr] measurements
% Tm ... [ntxnr] observed event-times
% ind_y ... [1xNY] indexing of measurements
% ind_t ... [1xNE] indexing of events
% F_diff ... [integer] number of derivatives of the hessian matrix
% b_diff ... [integer] number of derivatives of the estimated random effects
% s ... [integer] experimental index
% fms ...[boolean] flag indicating whe
%
% Outputs:
% ========
% estimated random effects bhat and hessian G and derivatives wrt to 
% beta and delta 
% the ordering of the output will depend on F_diff and b_diff
% see the file for details
%
% 2015/04/14 Fabian Froehlich



function varargout = optimize_SCTL_si(Model,Data,bhat_0,beta,delta,type_D,t,Ym,Tm,ind_y,ind_t,F_diff,b_diff,s,fms)
options_fmincon = optimset('algorithm','trust-region-reflective',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',100,... %'display','iter',...
    'TolFun',0,...
    'TolX',1e-10,...
    'PrecondBandWidth',Inf,...
    'Hessian','user-supplied');

if(fms)
    [bhat,OBJ,~,~,~,~,~] = fmincon(...
        @(b) objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),...
        bhat_0,[],[],[],[],-5*ones(length(bhat_0),1),5*ones(length(bhat_0),1),[],options_fmincon);
    rng(0);
    bhat_0_lhc = 10*lhsdesign(10,length(bhat_0),'smooth','off')' - 5;
    
    for j = 1:10
        try 
        [bhatp,OBJp,~,~,~,~,~] = fmincon(...
            @(b) objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),...
            bhat_0_lhc(:,j),[],[],[],[],-5*ones(length(bhat_0),1),5*ones(length(bhat_0),1),[],options_fmincon);
        if(OBJp<OBJ)
            bhat = bhatp;
        end
        catch
        end
    end
else
    [bhat,~,~,~,~,~,~] = fmincon(...
    @(b) objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),...
    bhat_0,[],[],[],[],-5*ones(length(bhat_0),1),5*ones(length(bhat_0),1),[],options_fmincon);
end


switch(b_diff)
    case 0
        % Higher order derivatives objective function
        [~,~,G] = objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
        varargout{1} = bhat;
        varargout{2} = G;
    case 1
        switch(F_diff)
            case 3
                % Higher order derivatives objective function
                [~,~,G,ddJdbdb,...
                    ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                    dGdb,pdGpdbeta,pdGpddelta] = ...
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
                % Gradient of optimal point
                [bhat,dbhatdbeta,dbhatddelta] = bhat_SCTL_si(bhat,ddJdbdb,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta);
                
                varargout{1} = bhat;
                
                varargout{2} = dbhatdbeta;
                varargout{3} = dbhatddelta;
                
                varargout{4} = G;
                
                varargout{5} = dGdb;
                varargout{6} = pdGpdbeta;
                varargout{7} = pdGpddelta;
            case 2
                % Higher order derivatives objective function
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,1,2)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,2,4)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(beta,@(beta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,2,5)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(delta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,2,6)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,3,10)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(beta,@(beta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,3,11)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(delta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,3,12)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(beta,@(beta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,4,13)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(delta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,5,14)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(delta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-5,4,15)
                [~,~,G,ddJdbdb,...
                    ddJdbdbeta,ddJdbddelta,~,~,~] = ...
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
                % Gradient of optimal point
                [bhat,dbhatdbeta,dbhatddelta] = bhat_SCTL_si(bhat,ddJdbdb,ddJdbdbeta,ddJdbddelta,[],[],[]);
                varargout{1} = bhat;
                varargout{2} = dbhatdbeta;
                varargout{3} = dbhatddelta;
                varargout{4} = G;
        end
    case 2
        switch(F_diff)
            case 4
                % Higher order derivatives objective function
                [~,~,G,ddJdbdb,...
                    ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                    dGdb,pdGpdbeta,pdGpddelta,...
                    dddJdbdbetadbeta,dddJdbddeltaddelta,dddJdbdbetaddelta,...
                    ddGdbdb,pddGdbpdbeta,pdpdGpdbetapdbeta,pddGdbpddelta,pdpdGpddeltapddelta,pdpdGpdbetapddelta] = ...
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
                % Gradient + Hessian of optimal point
                [bhat,dbhatdbeta,dbhatddelta,ddbhatdbetadbeta,ddbhatdbetaddelta,ddbhatddeltaddelta] = bhat_SCTL_si(bhat,ddJdbdb,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,dGdb,pdGpdbeta,pdGpddelta,dddJdbdbetadbeta,dddJdbddeltaddelta,dddJdbdbetaddelta);
                
                
                varargout{1} = bhat;
                
                varargout{2} = dbhatdbeta;
                varargout{3} = dbhatddelta;
                
                varargout{4} = ddbhatdbetadbeta;
                varargout{5} = ddbhatdbetaddelta;
                varargout{6} = ddbhatddeltaddelta;
                
                varargout{7} = G;
                
                varargout{8} = dGdb;
                varargout{9} = pdGpdbeta;
                varargout{10} = pdGpddelta;
                
                varargout{11} = ddGdbdb;
                varargout{12} = pddGdbpdbeta;
                varargout{13} = pdpdGpdbetapdbeta;
                varargout{14} = pddGdbpddelta;
                varargout{15} = pdpdGpddeltapddelta;
                varargout{16} = pdpdGpdbetapddelta;
                
            case 3
                % Higher order derivatives objective function
                [~,~,G,ddJdbdb,...
                    ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                    dGdb,pdGpdbeta,pdGpddelta,...
                    dddJdbdbetadbeta,dddJdbddeltaddelta,dddJdbdbetaddelta] = ...
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
                % Gradient + Hessian of estimated random effects
                [bhat,dbhatdbeta,dbhatddelta,ddbhatdbetadbeta,ddbhatdbetaddelta,ddbhatddeltaddelta] = bhat_SCTL_si(bhat,ddJdbdb,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,dGdb,pdGpdbeta,pdGpddelta,dddJdbdbetadbeta,dddJdbddeltaddelta,dddJdbdbetaddelta);
                
                varargout{1} = bhat;
                
                varargout{2} = dbhatdbeta;
                varargout{3} = dbhatddelta;
                
                varargout{4} = ddbhatdbetadbeta;
                varargout{5} = ddbhatdbetaddelta;
                varargout{6} = ddbhatddeltaddelta;
                
                varargout{7} = G;
                
                varargout{8} = dGdb;
                varargout{9} = pdGpdbeta;
                varargout{10} = pdGpddelta;
        end
        
        
end
end