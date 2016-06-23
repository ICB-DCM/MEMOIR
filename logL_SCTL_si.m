function [ logL_D, logL_T, logL_b, logL_I, bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i)
%LOGL_SCTL_SI Summary of this function goes here
%   Detailed explanation goes here

if(isfield(P{s},'SCTL'))
    if(isfield(P{s}.SCTL,'dbdxi'))
        bhat_si0 = P{s}.SCTL.bhat(:,i) + P{s}.SCTL.dbdxi(:,:,i)*(xi-P{s}.xi);
    else
        bhat_si0 = P{s}.SCTL.bhat(:,i);
    end
else
    bhat_si0 = zeros(length(model.ind_b),1);
end

% compute bhat (optimum of the single cell likelihood
% (objective_SCTL_s1) with respect to b and the respective
% derivatives with respect to beta and delta
%
% testing:
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) getBhat( beta, delta, bhat_si0, model, data, s, i, options, P),1e-5,'val','dbeta')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) getBhat( beta, delta, bhat_si0, model, data, s, i, options, P),1e-5,'val','ddelta')

beta = model.beta(xi);
delta = model.delta(xi);
if(options.nderiv > 0)
    dbetadxi = model.dbetadxi(xi);
    ddeltadxi = model.ddeltadxi(xi);
    ddbetadxidxi = model.ddbetadxidxi(xi);
    dddeltadxidxi = model.dddeltadxidxi(xi);
end

[B,G] = getBhat( beta, delta, bhat_si0, model, data, s, i, options, P);

% Construct single-cell parameter
phi_si = model.phi( beta, B.val);
Ym_si = data.SCTL.Y(:,:,i);
Tm_si = data.SCTL.T(:,:,i);
bhat_si = B.val;
t_s = data.SCTL.time(data.SCTL.ind_y(:,i));
ind_y = data.SCTL.ind_y(:,i);
ind_t = data.SCTL.ind_t(:,i);


% Simulate model and compute derivatives
if(options.nderiv == 0)
    [Y_si,T_si,R_si] = simulate_trajectory(t_s,phi_si,model,data.condition,s,ind_t,ind_y);
elseif(and(options.nderiv == 1,options.integration == 0))
    [Y_si,T_si,R_si,dY_sidphi,dT_sidphi,dR_sidphi] = simulate_trajectory(t_s,phi_si,model,data.condition,s,ind_t,ind_y);
else
    [Y_si,T_si,R_si,dY_sidphi,dT_sidphi,dR_sidphi,ddY_sidphidphi,ddT_sidphidphi,ddR_sidphidphi] = simulate_trajectory(t_s,phi_si,model,data.condition,s,ind_t,ind_y);
end

% Construct sigma
if(options.nderiv<1)
    [Sigma_noise_si] = build_sigma_noise(phi_si,Ym_si,s,model,ind_y);
    [Sigma_time_si] = build_sigma_time(phi_si,Tm_si,s,model,ind_t);
elseif(options.nderiv<2)
    [Sigma_noise_si,dSigma_noisedphi] = build_sigma_noise(phi_si,Ym_si,s,model,ind_y);
    [Sigma_time_si,dSigma_timedphi] = build_sigma_time(phi_si,Tm_si,s,model,ind_t);
else
    [Sigma_noise_si,dSigma_noisedphi,ddSigma_noisedphidphi] = build_sigma_noise(phi_si,Ym_si,s,model,ind_y);
    [Sigma_time_si,dSigma_timedphi,ddSigma_timedphidphi] = build_sigma_time(phi_si,Tm_si,s,model,ind_t);
end

%% Evaluation of likelihood and likelihood gradient

% this is part accounts for the noise model
% J_D = log(p(Y(b,beta)|D))
switch(model.noise_model)
    case 'normal'
        noisedist = @normal_noise;
    case 'lognormal'
        noisedist = @lognormal_noise;
    case 'tdist'
        noisedist = @tdist_noise;
end

switch(options.nderiv)
    case 0
        J_D = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
    case 1
        [J_D,dJ_DdY,dJ_DdSigma] = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
    case 2
        [J_D,dJ_DdY,dJ_DdSigma,ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
end

% this is part accounts for the event model
% J_D = log(p(Y(b,beta)|D))
if(~isfield(model,'time_model'))
    model.time_model = 'normal';
end
switch(model.time_model)
    case 'normal'
        switch(options.nderiv)
            case 0
                J_T = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
            case 1
                [J_T,dJ_TdT,dJ_TdR,dJ_TdSigma] = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
            case 2
                [J_T,dJ_TdT,dJ_TdR,dJ_TdSigma,ddJ_TdTdT,ddJ_TdTdR,ddJ_TdRdR,ddJ_TdTdSigma,ddJ_TdRdSigma,ddJ_TdSigmadSigma] = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
        end
end

% this part accounts for the parameter model
% J_b = log(p(b_si|delta))
switch(model.parameter_model)
    case 'normal'
        switch(options.nderiv)
            case 0
                J_b = normal_param(bhat_si,delta,options.type_D);
            case 1
                [J_b,dJ_bdb,pdJ_bpddelta]= normal_param(bhat_si,delta,options.type_D);
            case 2
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(bhat_si,@(b) normal_param(b,delta,options.type_D),1e-4,1,2)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,options.type_D),1e-4,1,3)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(bhat_si,@(b) normal_param(b,delta,options.type_D),1e-4,2,4)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,options.type_D),1e-4,2,5)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,options.type_D),1e-4,3,6)
                [J_b,dJ_bdb,pdJ_bpddelta,ddJ_bdbdb,dpdJ_bdbpddelta,pdpdJ_bpddeltapddelta]= normal_param(bhat_si,delta,options.type_D);
        end
    case 'lognormal'
        switch(options.nderiv)
            case 0
                J_b = lognormal_param(bhat_si,delta,options.type_D);
            case 1
                [J_b,dJ_bdb,pdJ_bpddelta]= lognormal_param(bhat_si,delta,options.type_D);
            case 2
                [J_b,dJ_bdb,pdJ_bpddelta,ddJ_bdbdb,dpdJ_bdbpddelta,pdpdJ_bpddeltapddelta] = lognormal_param(bhat_si,delta,options.type_D);
        end
end

logL_D.val =  - J_D;
logL_T.val  = - J_T;
logL_b.val  = - J_b;
bhat.val = B.val;

if(options.integration)
    % laplace approximation
    logL_I.val = - 0.5*log(det(G.val));
end

if options.nderiv >= 1
    dbhat_sidbeta = B.dbeta;
    dbhat_siddelta = B.ddelta;
    % first order derivatives
    dphidb = model.dphidb(beta,bhat_si);
    pdphipdbeta  = model.dphidbeta(beta,bhat_si);
    
    dphidbeta = chainrule(dphidb,dbhat_sidbeta) + pdphipdbeta;
    dphiddelta = chainrule(dphidb,dbhat_siddelta);
    dphidxi = chainrule(dphidbeta,dbetadxi) + chainrule(dphiddelta,ddeltadxi);
    
    dJ_Ddphi = chainrule(dJ_DdY,dY_sidphi) + chainrule(dJ_DdSigma,dSigma_noisedphi) ;
    dJ_Ddxi = chainrule(dJ_Ddphi,dphidxi);
    
    dJ_Tdphi = chainrule(dJ_TdT,dT_sidphi) + chainrule(dJ_TdR,dR_sidphi) + chainrule(dJ_TdSigma,dSigma_timedphi) ;
    dJ_Tdxi = chainrule(dJ_Tdphi,dphidxi);
    
    dbdxi = chainrule(dbhat_sidbeta,dbetadxi) + chainrule(dbhat_siddelta,ddeltadxi);
    bhat.dxi = dbdxi;
    
    dJ_bdxi = chainrule(dJ_bdb,dbdxi) + chainrule(pdJ_bpddelta,ddeltadxi);
    
    logL_D.dxi = - transpose(dJ_Ddxi);
    logL_T.dxi = - transpose(dJ_Tdxi);
    logL_b.dxi = - transpose(dJ_bdxi);
    
    if(options.integration)
        % laplace approximation
        invG = pinv(G.val);
        
        % take care when nelem(b) == 1 ... (ones(1,1,1) ==
        % ones(1,1) so G.db will be missing one dimension!)
        if(numel(bhat_si)==1)
            G.dbeta = pdGpdbeta + permute(G.db*dbhat_sidbeta,[3,1,2]);
            dGddelta = pdGpddelta + permute(G.db*dbhat_siddelta,[3,1,2]);
            dGdxi = chainrule(G.dbeta,dbetadxi) + chainrule(dGddelta,ddeltadxi);
        else
            G.dbeta = G.dbeta + chainrule(G.db,dbhat_sidbeta);
            dGddelta = G.ddelta + chainrule(G.db,dbhat_siddelta);
            dGdxi = chainrule(G.dbeta,dbetadxi) + chainrule(dGddelta,ddeltadxi);
        end
        
        logL_I.dxi = - 0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,invG,permute(dGdxi,[4,1,2,3])),2)),eye(length(bhat_si))),1),2)); % 1/2*Tr(invG*dG)
    end
    
    if options.nderiv >= 2
        % second order derivatives
        
        ddbhat_sidbetadbeta = B.dbetadbeta;
        ddbhat_sidbetaddelta = B.dbetaddelta;
        ddbhat_siddeltaddelta = B.ddeltaddelta;
        
        ddphidbdb = model.ddphidbdb(beta,bhat_si);
        ddphidbdbeta = model.ddphidbdbeta(beta,bhat_si);
        
        ddphidbetadbeta = chainrule(dphidb,ddbhat_sidbetadbeta) + chainrule_ddxdydy_dydz(ddphidbdb,dbhat_sidbeta) ...
            + permute(chainrule(permute(ddphidbdbeta,[1,3,2]),dbhat_sidbeta),[1,3,2]);
        ddphidbetaddelta = chainrule(dphidb,ddbhat_sidbetaddelta) + chainrule_ddxdydy_dydz_dydv(ddphidbdb,dbhat_sidbeta,dbhat_siddelta);
        ddphiddeltaddelta = chainrule(dphidb,ddbhat_siddeltaddelta) + chainrule_ddxdydy_dydz(ddphidbdb,dbhat_siddelta);
        
        if(numel(bhat_si)==1) % we have to do this manually here since 3rd order tensors with trailing 1 dimensional orders are not possible in matlab ...
            ddphidxidxi = chainrule(dphidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddphidbetadbeta,dbetadxi) ...
                + chainrule(dphiddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddphiddeltaddelta,ddeltadxi) ...
                + permute(sum(bsxfun(@times,bsxfun(@times,ddphidbetaddelta,permute(dbetadxi,[3,1,2])),permute(ddeltadxi,[3,4,5,2,1])),2),[1,3,4,2]) ...
                + permute(sum(bsxfun(@times,bsxfun(@times,ddphidbetaddelta,permute(dbetadxi,[3,1,2])),permute(ddeltadxi,[3,4,5,2,1])),2),[1,4,3,2]);
        else
            ddphidxidxi = chainrule(dphidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddphidbetadbeta,dbetadxi) ...
                + chainrule(dphiddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddphiddeltaddelta,ddeltadxi) ...
                + chainrule_ddxdydy_dydz_dydv(ddphidbetaddelta,dbetadxi,ddeltadxi) + permute(chainrule_ddxdydy_dydz_dydv(ddphidbetaddelta,dbetadxi,ddeltadxi),[1,3,2]);
        end
        
        ddJ_DdphidY = bsxfun(@times,ddJ_DdYdY,permute(dY_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_DdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
        
        ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dY_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
        
        ddJ_Ddphidphi = chainrule(dJ_DdY,ddY_sidphidphi) ...
            + squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(dY_sidphi,[3,1,4,2])) ...
            + bsxfun(@times,ddJ_DdphidSigma,permute(dSigma_noisedphi,[3,1,4,2])),2));
        
        ddJ_Ddxidxi = chainrule(dJ_Ddphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Ddphidphi,dphidxi);
        
        ddJ_TdphidT = bsxfun(@times,ddJ_TdTdT,permute(dT_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_TdTdR,permute(dR_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_TdTdSigma,permute(dSigma_timedphi,[3,1,2]));
        
        ddJ_TdphidR = bsxfun(@times,ddJ_TdTdR,permute(dR_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_TdRdR,permute(dR_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_TdRdSigma,permute(dSigma_timedphi,[3,1,2]));
        
        ddJ_TdphidSigma = bsxfun(@times,ddJ_TdTdSigma,permute(dT_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_TdRdSigma,permute(dR_sidphi,[3,1,2])) + ...
            bsxfun(@times,ddJ_TdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
        
        ddJ_Tdphidphi = chainrule(dJ_TdT,ddT_sidphidphi) ...
            + chainrule(dJ_TdR,ddR_sidphidphi) ...
            + squeeze(sum(bsxfun(@times,ddJ_TdphidT,permute(dT_sidphi,[3,1,4,2])) ...
            + bsxfun(@times,ddJ_TdphidR,permute(dR_sidphi,[3,1,4,2])) ...
            + bsxfun(@times,ddJ_TdphidSigma,permute(dSigma_timedphi,[3,1,4,2])),2));
        
        ddJ_Tdxidxi = chainrule(dJ_Tdphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Tdphidphi,dphidxi);
        
        ddbdxidxi = chainrule(dbhat_sidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddbhat_sidbetadbeta,dbetadxi) ...
            + chainrule(dbhat_siddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddbhat_siddeltaddelta,ddeltadxi) ...
            + chainrule_ddxdydy_dydz_dydv(ddbhat_sidbetaddelta,dbetadxi,ddeltadxi) ...
            + chainrule_ddxdydy_dydz_dydv(permute(ddbhat_sidbetaddelta,[1,3,2]),ddeltadxi,dbetadxi);
        
        bhat.dxidxi = ddbdxidxi;
        
        ddJ_bdxidxi = chainrule(dJ_bdb,ddbdxidxi) + chainrule_ddxdydy_dydz(ddJ_bdbdb,dbdxi) ...
            + chainrule(pdJ_bpddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(pdpdJ_bpddeltapddelta,ddeltadxi) ...
            + chainrule_ddxdydy_dydz_dydv(dpdJ_bdbpddelta,dbdxi,ddeltadxi) ...
            + chainrule_ddxdydy_dydz_dydv(permute(dpdJ_bdbpddelta,[2,1]),ddeltadxi,dbdxi);
        
        logL_D.dxidxi = - ddJ_Ddxidxi;
        logL_T.dxidxi = - ddJ_Tdxidxi;
        logL_b.dxidxi = - ddJ_bdxidxi;
        
        if(options.integration)
            % laplace approximation
            invG = pinv(G.val);
            
            ddGdbetadbeta = G.dbetadbeta + 2*chainrule(permute(G.dbdbeta,[1,2,4,3]),dbhat_sidbeta) ...
                + chainrule_ddxdydy_dydz(G.dbdb,dbhat_sidbeta) + chainrule(G.db,ddbhat_sidbetadbeta);
            
            ddGddeltaddelta = G.ddeltaddelta + 2*chainrule(permute(G.dbddelta,[1,2,4,3]),dbhat_siddelta) ...
                + chainrule_ddxdydy_dydz(G.dbdb,dbhat_siddelta) + chainrule(G.db,ddbhat_siddeltaddelta);
            
            ddGdbetaddelta = G.dbetaddelta + chainrule(permute(G.dbdbeta,[1,2,4,3]),dbhat_siddelta) ...
                + permute(chainrule(permute(G.dbddelta,[1,2,4,3]),dbhat_sidbeta),[1,2,4,3]) ...
                + chainrule_ddxdydy_dydz_dydv(G.dbdb,dbhat_sidbeta,dbhat_siddelta) ...
                + chainrule(G.db,ddbhat_sidbetaddelta);
            
            ddGdxidxi = chainrule_ddxdydy_dydz(ddGdbetadbeta,dbetadxi) + chainrule(G.dbeta,ddbetadxidxi) ...
                + chainrule_ddxdydy_dydz(ddGddeltaddelta,ddeltadxi) + chainrule(dGddelta,dddeltadxidxi) ...
                + 2*chainrule_ddxdydy_dydz_dydv(ddGdbetaddelta,dbetadxi,ddeltadxi);
            
            dinvGdxi = squeeze(sum(bsxfun(@times,invG,permute(squeeze(sum(bsxfun(@times,permute(dGdxi,[4,1,2,3]),invG),2)),[4,1,2,3])),2));
            
            logL_I.dxidxi = - 0.5*squeeze(sum(sum(squeeze(bsxfun(@times,sum(bsxfun(@times,permute(dinvGdxi,[1,2,4,3]),permute(dGdxi,[4,1,2,5,3])),2)+sum(bsxfun(@times,invG,permute(ddGdxidxi,[5,1,2,3,4])),2),permute(eye(length(bhat_si)),[1,3,2]))),1),2)); % 1/2*Tr(dinvG*dg + invG*ddG)
        end
        
    end
    
end

Sim.SCTL_Y = reshape(Y_si,[sum(ind_y),size(data.SCTL.Y,2)]);
Sim.SCTL_T = reshape(T_si(ind_t),[sum(ind_t),size(data.SCTL.T,2)]);
Sim.SCTL_R = reshape(R_si(ind_t),[sum(ind_t),size(data.SCTL.T,2)]);


end

