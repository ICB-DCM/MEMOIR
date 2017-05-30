% optimize_SCTL_si estimates the random effects parameters and computes
% derivatives wrt beta and delta moreover some derivatives of the hessian
% of the respective objective function are also returned as they are needed
% for certain computations later on and computational complexity is rather
% high
%
% USAGE:
% ======
% [...] = bhat_SCTL_si(model,data,bhat_0,beta,delta,type_D,t,Ym,Tm,ind_y,ind_t,F_diff,b_diff,s,fms)
%
% INPUTS:
% =======
% model ... [struct] model definition
% data ... [struct] data definition
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



function [B,FIM,J,Sim] = optimize_SCTL_si(model,data,bhat_0,beta,delta,F_diff,b_diff,s,i,options,P_old)
options_fmincon = optimset('algorithm','trust-region-reflective',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',1000,...
    'TolFun',0,...
    'TolX',1e-10,...
    'PrecondBandWidth',Inf,...
    'Hessian','user-supplied');

    % do multistart every few iterations
    fms = (mod(P_old{s}.n_store,options.ms_iter)==0);

if(fms)
    try
    [bhat,OBJ,~,~,~,~,~] = fmincon(...
        @(b) objective_SCTL_s1(model,data,beta,b,delta,s,i,options,1),...
        bhat_0,[],[],[],[],-10*ones(length(bhat_0),1),10*ones(length(bhat_0),1),[],options_fmincon);
    catch err
        OBJ = inf;
    end
    rng(0);
    N_MS = 20;
    bhat_0_lhc = [bhat_0,6*lhsdesign(N_MS,length(bhat_0),'smooth','off')' - 3];
    OBJs = zeros(N_MS+1,1);
    OBJs(1) = OBJ;

    for j = 1:N_MS
        try 
        [bhatp,OBJp,~,~,~,~,~] = fmincon(...
            @(b) objective_SCTL_s1(model,data,beta,b,delta,s,i,options,1),...
            bhat_0_lhc(:,j),[],[],[],[],-10*ones(length(bhat_0),1),10*ones(length(bhat_0),1),[],options_fmincon);
        OBJs(1+j) = OBJp;
        if(OBJp<OBJ)
            bhat = bhatp;
            OBJ = OBJp;
        end
        catch
        end
    end
else
    [bhat,OBJ,~,~,~,~,~] = fmincon(...
    @(b) objective_SCTL_s1(model,data,beta,b,delta,s,i,options,1),...
    bhat_0,[],[],[],[],-10*ones(length(bhat_0),1),10*ones(length(bhat_0),1),[],options_fmincon);
end

% testing:
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(bhat,@(bhat) objSCTL_J(model,data,beta,bhat,delta,s,i,options,1),1e-3,'val','db',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,1),1e-3,'val','dbeta',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,1),1e-3,'val','ddelta',true)
%
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(bhat,@(bhat) objSCTL_J(model,data,beta,bhat,delta,s,i,options,2),1e-3,'db','dbdb',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,2),1e-3,'db','dbdbeta',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,2),1e-3,'db','dbddelta',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,2),1e-3,'dbeta','dbetadbeta',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,2),1e-3,'ddelta','ddeltaddelta',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) objSCTL_J(model,data,beta,bhat,delta,s,i,options,2),1e-3,'dbeta','dbetaddelta',true)
%
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(bhat,@(bhat) objSCTL_FIM(model,data,beta,bhat,delta,s,i,options,3),1e-4,'val','db',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) objSCTL_FIM(model,data,beta,bhat,delta,s,i,options,3),1e-3,'val','dbeta',true)
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) objSCTL_FIM(model,data,beta,bhat,delta,s,i,options,3),1e-6,'val','ddelta',true)

assert(OBJ<inf);

[~,~,~,J,FIM,Sim] = objective_SCTL_s1(model,data,beta,bhat,delta,s,i,options,F_diff);
B = bhat_SCTL_si(bhat,FIM,J,b_diff);

end

%% functions to test gradient

function J = objSCTL_J(model,data,beta,bhat,delta,s,i,options,F_diff)
    [~,~,~,J] = objective_SCTL_s1(model,data,beta,bhat,delta,s,i,options,F_diff);
end

function FIM = objSCTL_FIM(model,data,beta,bhat,delta,s,i,options,F_diff)
    [~,~,~,~,FIM] = objective_SCTL_s1(model,data,beta,bhat,delta,s,i,options,F_diff);
end