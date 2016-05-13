function [P,logL_sc,dlogL_scdxi,ddlogL_scdxidxi] = logL_SCTL(xi, model, data, s, options, P)
    
    persistent fp
    persistent fl
    
    options.nderiv = nargout-2;
    
    %% Construct fixed effects and covariance matrix
    beta = model.beta(xi);
    delta = model.delta(xi);
    
    [D,~,~,~,~,~] = xi2D(delta,options.type_D);
    % debugging:
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,1,3)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,3,5)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,2,4)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,4,6)
    
    % Initialization measurements
    Sim_SCTL.Y = NaN(size(data.SCTL.Y));
    % events
    if(~isfield(data.SCTL,'T'))
        data.SCTL.T = zeros(0,1,size(data.SCTL.Y,3));
    end
    Sim_SCTL.T = NaN(size(data.SCTL.T));
    Sim_SCTL.R = NaN(size(data.SCTL.T));
    
    % set default scaling
    if(~isfield(model,'SCTLscale'))
        model.SCTLscale = 1;
    end
    
    % Loop: Indiviudal cells
    
    dbetadxi = model.dbetadxi(xi);
    ddeltadxi = model.ddeltadxi(xi);
    ddbetadxidxi = model.ddbetadxidxi(xi);
    dddeltadxidxi = model.dddeltadxidxi(xi);
    
    logLi_D = zeros(1,size(data.SCTL.Y,3));
    logLi_T = zeros(1,size(data.SCTL.Y,3));
    logLi_b = zeros(1,size(data.SCTL.Y,3));
    logLi_I = zeros(1,size(data.SCTL.Y,3));
    if options.nderiv >= 1
        dlogLi_Ddxi = zeros(length(xi),size(data.SCTL.Y,3));
        dlogLi_Tdxi = zeros(length(xi),size(data.SCTL.Y,3));
        dlogLi_bdxi = zeros(length(xi),size(data.SCTL.Y,3));
        dlogLi_Idxi = zeros(length(xi),size(data.SCTL.Y,3));
        if options.nderiv > 2
            ddlogLi_Ddxidxi = zeros(length(xi),length(xi),size(data.SCTL.Y,3));
            ddlogLi_Tdxidxi = zeros(length(xi),length(xi),size(data.SCTL.Y,3));
            ddlogLi_bdxidxi = zeros(length(xi),length(xi),size(data.SCTL.Y,3));
            ddlogLi_Idxidxi = zeros(length(xi),length(xi),size(data.SCTL.Y,3));
        end
    end

    tmp = arrayfun(@(x) any(~isnan(data.SCTL.Y(:,:,x)),2),1:size(data.SCTL.Y,3),'UniformOutput',false);
    data.SCTL.ind_y = [tmp{:}];
    tmp = arrayfun(@(x) any(~isnan(data.SCTL.T(:,:,x)),2),1:size(data.SCTL.Y,3),'UniformOutput',false);
    data.SCTL.ind_t = [tmp{:}];
    
    beta = model.beta(xi);
    delta = model.delta(xi);
    
    for i = 1:size(data.SCTL.Y,3)

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
        
        logLi_D(1,i) =  - J_D;
        logLi_T(1,i)  = - J_T;
        logLi_b(1,i)  = - J_b;
        bhat(:,i) = B.val;
        
        if(options.integration)
            % laplace approximation
            logLi_I(1,i) = - 0.5*log(det(G.val));
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
            dbhatdxi(:,:,i) = dbdxi;
            
            dJ_bdxi = chainrule(dJ_bdb,dbdxi) + chainrule(pdJ_bpddelta,ddeltadxi);
            
            dlogLi_Ddxi(:,i) = - transpose(dJ_Ddxi);
            dlogLi_Tdxi(:,i) = - transpose(dJ_Tdxi);
            dlogLi_bdxi(:,i) = - transpose(dJ_bdxi);
            
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
                
                dlogLi_Idxi(:,i) = - 0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,invG,permute(dGdxi,[4,1,2,3])),2)),eye(length(bhat_si))),1),2)); % 1/2*Tr(invG*dG)
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
                
                ddbhatdxidxi(:,:,:,i) = ddbdxidxi;
                
                ddJ_bdxidxi = chainrule(dJ_bdb,ddbdxidxi) + chainrule_ddxdydy_dydz(ddJ_bdbdb,dbdxi) ...
                    + chainrule(pdJ_bpddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(pdpdJ_bpddeltapddelta,ddeltadxi) ...
                    + chainrule_ddxdydy_dydz_dydv(dpdJ_bdbpddelta,dbdxi,ddeltadxi) ...
                    + chainrule_ddxdydy_dydz_dydv(permute(dpdJ_bdbpddelta,[2,1]),ddeltadxi,dbdxi);
                
                ddlogLi_Ddxidxi(:,:,i) = - ddJ_Ddxidxi;
                ddlogLi_Tdxidxi(:,:,i) = - ddJ_Tdxidxi;
                ddlogLi_bdxidxi(:,:,i) = - ddJ_bdxidxi;
                
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
                    
                    ddlogLi_Idxidxi(:,:,i) = - 0.5*squeeze(sum(sum(squeeze(bsxfun(@times,sum(bsxfun(@times,permute(dinvGdxi,[1,2,4,3]),permute(dGdxi,[4,1,2,5,3])),2)+sum(bsxfun(@times,invG,permute(ddGdxidxi,[5,1,2,3,4])),2),permute(eye(length(bhat_si)),[1,3,2]))),1),2)); % 1/2*Tr(dinvG*dg + invG*ddG)
                end
                
            end
            
        end
        
        Sim_SCTL.Y(ind_y,:,i) = reshape(Y_si,[sum(ind_y),size(Sim_SCTL.Y,2)]);
        Sim_SCTL.T(ind_t,:,i) = reshape(T_si,[sum(ind_t),size(Sim_SCTL.T,2)]);
        Sim_SCTL.R(ind_t,:,i) = reshape(R_si,[sum(ind_t),size(Sim_SCTL.R,2)]);
    end
    
    P{s}.SCTL.bhat = bhat;
    if(options.nderiv>0)
      P{s}.SCTL.dbhatdxi = dbhatdxi;  
    end
    
    %% Visulization
    if options.plot
        
        % Visualisation of single cell parameters
        if(isempty(fp))
            if(isfield(model,'title'))
                if(ischar(model.title))
                    fp(s) = figure('Name',model.title);
                else
                    fp(s) = figure;
                end
            else
                fp(s) = figure;
            end
        else
            if(length(fp)<s)
                if(isfield(model,'title'))
                    if(ischar(model.title))
                        fp(s) = figure('Name',model.title);
                    else
                        fp(s) = figure;
                    end
                else
                    fp(s) = figure;
                end
            elseif(isempty(fp(s)))
                if(isfield(model,'title'))
                    if(ischar(model.title))
                        fp(s) = figure('Name',model.title);
                    else
                        fp(s) = figure;
                    end
                else
                    fp(s) = figure;
                end
            end
        end
        figure(fp(s))
        clf
        b_s = P{s}.SCTL.bhat;
        n_b = size(b_s,1);
        
        for j = 1:n_b
            subplot(ceil((n_b+1)/4),4,j+1)
            xx = linspace(-5*sqrt(D(j,j)),5*sqrt(D(j,j)),100);
            %nhist(P{s}.SCTL.bhat(j,:),'pdf','noerror');
            hold on
            plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','LineWidth',2)
            ecdf = zeros(length(xx),1);
            for k = 1:length(xx)
                ecdf(k) = sum(b_s(j,:)<xx(k))/length(b_s(j,:));
            end
            plot(xx,ecdf,'--r','LineWidth',2)
            
            
            if(j==1)
                
            end
            xlim([-5*sqrt(D(j,j)),5*sqrt(D(j,j))])
            ylim([0,1.1])
            %xlabel(char(model.sym.b(model.ind_b(j))));
            ylabel('cdf')
            box on
        end
        subplot(ceil(n_b+1/4),4,1,'Visible','off')
        hold on
        plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','Visible','off')
        plot(xx,ecdf,'--r','LineWidth',2,'Visible','off')
        
        
        legend('cdf of single cell Parameters','cdf of populaton Parameters')
        
        % Visualisation of likelihood contribution
        if(isempty(fl))
            if(isfield(model,'title'))
                if(ischar(model.title))
                    fl(s) = figure('Name',model.title);
                else
                    fl(s) = figure;
                end
            else
                fl(s) = figure;
            end
        else
            if(length(fl)<s)
                if(isfield(model,'title'))
                    if(ischar(model.title))
                        fl(s) = figure('Name',model.title);
                    else
                        fl(s) = figure;
                    end
                else
                    fl(s) = figure;
                end
            elseif(isempty(fl(s)))
                if(isfield(model,'title'))
                    if(ischar(model.title))
                        fl(s) = figure('Name',model.title);
                    else
                        fl(s) = figure;
                    end
                else
                    fl(s) = figure;
                end
            end
        end
        figure(fl(s))
        clf
        if(options.integration)
            bar([logLi_D;logLi_T;logLi_b;logLi_I],'stacked')
            set(gca,'XTickLabel',{'data','event','par','int'})
        else
            bar([logLi_D;logLi_T;logLi_b],'stacked')
            set(gca,'XTickLabel',{'data','event','par'})
        end
        ylabel('log-likelihood')
        title('likelihood contribution')
        
        
        
        % Visualisation of data and fit
        model.plot(data,Sim_SCTL,s);
    end
    
    %% Summation
    logL_sc = logLi_D + logLi_b;
    if(options.events)
        logL_sc = logL_sc + logLi_T;
    end
    if(options.integration)
        logL_sc = logL_sc + logLi_I;
    end
    if(options.nderiv >= 1)
        dlogL_scdxi = dlogLi_Ddxi + dlogLi_bdxi;
        if(options.events)
            dlogL_scdxi = dlogL_scdxi + dlogLi_Tdxi;
        end
        if(options.integration)
            dlogL_scdxi = dlogL_scdxi + dlogLi_Idxi;
        end
        
        if(options.nderiv >= 2)
            ddlogL_scdxidxi = ddlogLi_Ddxidxi + ddlogLi_bdxidxi;
            if(options.events)
                ddlogL_scdxidxi = ddlogL_scdxidxi + ddlogLi_Tdxidxi;
            end
            if(options.integration)
                ddlogL_scdxidxi = ddlogL_scdxidxi + ddlogLi_Idxidxi;
            end
        end
    end
    
    
end
