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



function [B,J,FIM] = optimize_SCTL_si(Model,Data,bhat_0,beta,delta,type_D,t,Ym,Tm,ind_y,ind_t,F_diff,b_diff,s,fms)
options_fmincon = optimset('algorithm','trust-region-reflective',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',100,... % 'display','iter',...
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

[J,FIM] = objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s,F_diff);
B = bhat_SCTL_si(bhat,FIM,J,b_diff);

end