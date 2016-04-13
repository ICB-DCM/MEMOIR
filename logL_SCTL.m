function [P,logL_SCTL,dlogL_SCTLdxi,ddlogL_SCTLdxi2] = logL_SCTL(xi, Model, Data, s, options, P, P_old)
    
    %% Construct fixed effects and covariance matrix
    beta = Model.exp{s}.beta(xi);
    delta = Model.exp{s}.delta(xi);
    
    [D,~,~,~,~,~] = xi2D(delta,type_D);
    % debugging:
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,1,3)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,3,5)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,2,4)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,4,6)
    
    % Initialization measurements
    Sim_SCTL.Y = NaN(size(Data{s}.SCTL.Y));
    % events
    if(~isfield(Data{s}.SCTL,'T'))
        Data{s}.SCTL.T = zeros(0,1,size(Data{s}.SCTL.Y,3));
    end
    Sim_SCTL.T = NaN(size(Data{s}.SCTL.T));
    Sim_SCTL.R = NaN(size(Data{s}.SCTL.T));
    
    % set default scaling
    if(~isfield(Model,'SCTLscale'))
        Model.SCTLscale = 1;
    end
    
    % Loop: Indiviudal cells
    
    dbetadxi = Model.exp{s}.dbetadxi(xi);
    ddeltadxi = Model.exp{s}.ddeltadxi(xi);
    ddbetadxidxi = Model.exp{s}.ddbetadxidxi(xi);
    dddeltadxidxi = Model.exp{s}.dddeltadxidxi(xi);
    
    logLi_D = zeros(1,size(Data{s}.SCTL.Y,3));
    logLi_T = zeros(1,size(Data{s}.SCTL.Y,3));
    logLi_b = zeros(1,size(Data{s}.SCTL.Y,3));
    logLi_I = zeros(1,size(Data{s}.SCTL.Y,3));
    if nderiv >= 1
        dlogLi_Ddxi = zeros(length(xi),size(Data{s}.SCTL.Y,3));
        dlogLi_Tdxi = zeros(length(xi),size(Data{s}.SCTL.Y,3));
        dlogLi_bdxi = zeros(length(xi),size(Data{s}.SCTL.Y,3));
        dlogLi_Idxi = zeros(length(xi),size(Data{s}.SCTL.Y,3));
        if nderiv > 2
            ddlogLi_Ddxidxi = zeros(length(xi),length(xi),size(Data{s}.SCTL.Y,3));
            ddlogLi_Tdxidxi = zeros(length(xi),length(xi),size(Data{s}.SCTL.Y,3));
            ddlogLi_bdxidxi = zeros(length(xi),length(xi),size(Data{s}.SCTL.Y,3));
            ddlogLi_Idxidxi = zeros(length(xi),length(xi),size(Data{s}.SCTL.Y,3));
        end
    end
    
    for i = 1:size(Data{s}.SCTL.Y,3)

        % Load single-cell data
        Ym_si = Data{s}.SCTL.Y(:,:,i);
        Tm_si = Data{s}.SCTL.T(:,:,i);
        

        

        
        
        
        % Construct single-cell parameter
        phi_si = Model.exp{s}.phi(beta,bhat_si);
        
        % Simulate model and compute derivatives
        if(nderiv == 1)
            [Y_si,T_si,R_si] = simulate_trajectory(t_s,phi_si,Model,Data{s}.condition,s,ind_t,ind_y);
        elseif(and(nderiv == 2,Model.integration == 0))
            [Y_si,T_si,R_si,dY_sidphi,dT_sidphi,dR_sidphi] = simulate_trajectory(t_s,phi_si,Model,Data{s}.condition,s,ind_t,ind_y);
        else
            [Y_si,T_si,R_si,dY_sidphi,dT_sidphi,dR_sidphi,ddY_sidphidphi,ddT_sidphidphi,ddR_sidphidphi] = simulate_trajectory(t_s,phi_si,Model,Data{s}.condition,s,ind_t,ind_y);
        end
        
        % Construct sigma
        if(nderiv<2)
            [Sigma_noise_si] = build_sigma_noise(phi_si,Ym_si,s,Model,ind_y);
            [Sigma_time_si] = build_sigma_time(phi_si,Tm_si,s,Model,ind_t);
        elseif(nderiv<3)
            [Sigma_noise_si,dSigma_noisedphi] = build_sigma_noise(phi_si,Ym_si,s,Model,ind_y);
            [Sigma_time_si,dSigma_timedphi] = build_sigma_time(phi_si,Tm_si,s,Model,ind_t);
        else
            [Sigma_noise_si,dSigma_noisedphi,ddSigma_noisedphidphi] = build_sigma_noise(phi_si,Ym_si,s,Model,ind_y);
            [Sigma_time_si,dSigma_timedphi,ddSigma_timedphidphi] = build_sigma_time(phi_si,Tm_si,s,Model,ind_t);
        end
        
        %% Evaluation of likelihood and likelihood gradient
        
        % this is part accounts for the noise model
        % J_D = log(p(Y(b,beta)|D))
        switch(Model.exp{s}.noise_model)
            case 'normal'
                noisedist = @normal_noise;
            case 'lognormal'
                noisedist = @lognormal_noise;
            case 'tdist'
                noisedist = @tdist_noise;
        end
        
        switch(nderiv)
            case 0
                J_D = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
            case 1
                J_D = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
            case 2
                [J_D,dJ_DdY,dJ_DdSigma] = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
            case 3
                [J_D,dJ_DdY,dJ_DdSigma,ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = noisedist(Y_si,Ym_si,Sigma_noise_si,ind_y);
        end
        
        % this is part accounts for the event model
        % J_D = log(p(Y(b,beta)|D))
        if(~isfield(Model.exp{s},'time_model'))
            Model.exp{s}.time_model = 'normal';
        end
        switch(Model.exp{s}.time_model)
            case 'normal'
                switch(nderiv)
                    case 0
                        J_T = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                    case 1
                        J_T = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                    case 2
                        [J_T,dJ_TdT,dJ_TdR,dJ_TdSigma] = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                    case 3
                        [J_T,dJ_TdT,dJ_TdR,dJ_TdSigma,ddJ_TdTdT,ddJ_TdTdR,ddJ_TdRdR,ddJ_TdTdSigma,ddJ_TdRdSigma,ddJ_TdSigmadSigma] = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                end
        end
        
        % this part accounts for the parameter model
        % J_b = log(p(b_si|delta))
        switch(Model.exp{s}.parameter_model)
            case 'normal'
                switch(nderiv)
                    case 0
                        J_b = normal_param(bhat_si,delta,type_D);
                    case 1
                        J_b = normal_param(bhat_si,delta,type_D);
                    case 2
                        [J_b,dJ_bdb,pdJ_bpddelta]= normal_param(bhat_si,delta,type_D);
                    case 3
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(bhat_si,@(b) normal_param(b,delta,type_D),1e-4,1,2)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,type_D),1e-4,1,3)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(bhat_si,@(b) normal_param(b,delta,type_D),1e-4,2,4)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,type_D),1e-4,2,5)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,type_D),1e-4,3,6)
                        [J_b,dJ_bdb,pdJ_bpddelta,ddJ_bdbdb,dpdJ_bdbpddelta,pdpdJ_bpddeltapddelta]= normal_param(bhat_si,delta,type_D);
                end
            case 'lognormal'
                switch(nderiv)
                    case 0
                        J_b = lognormal_param(bhat_si,delta,type_D);
                    case 1
                        J_b = lognormal_param(bhat_si,delta,type_D);
                    case 2
                        [J_b,dJ_bdb,pdJ_bpddelta]= lognormal_param(bhat_si,delta,type_D);
                    case 3
                        [J_b,dJ_bdb,pdJ_bpddelta,ddJ_bdbdb,dpdJ_bdbpddelta,pdpdJ_bpddeltapddelta] = lognormal_param(bhat_si,delta,type_D);
                end
        end
        
        logLi_D(1,i) =  - J_D;
        logLi_T(1,i)  = - J_T;
        logLi_b(1,i)  = - J_b;
        
        if(Model.integration)
            % laplace approximation
            logLi_I(1,i) = - 0.5*log(det(G));
        end
        
        if nderiv >= 2
            % first order derivatives
            dphidb = Model.exp{s}.dphidb(beta,bhat_si);
            pdphipdbeta  = Model.exp{s}.dphidbeta(beta,bhat_si);
            
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
            
            if(Model.integration)
                % laplace approximation
                invG = pinv(G);
                
                % take care when nelem(b) == 1 ... (ones(1,1,1) ==
                % ones(1,1) so dGdb will be missing one dimension!)
                if(numel(bhat_si)==1)
                    dGdbeta = pdGpdbeta + permute(dGdb*dbhat_sidbeta,[3,1,2]);
                    dGddelta = pdGpddelta + permute(dGdb*dbhat_siddelta,[3,1,2]);
                    dGdxi = chainrule(dGdbeta,dbetadxi) + chainrule(dGddelta,ddeltadxi);
                else
                    dGdbeta = pdGpdbeta + chainrule(dGdb,dbhat_sidbeta);
                    dGddelta = pdGpddelta + chainrule(dGdb,dbhat_siddelta);
                    dGdxi = chainrule(dGdbeta,dbetadxi) + chainrule(dGddelta,ddeltadxi);
                end
                
                dlogLi_Idxi(:,i) = - 0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,invG,permute(dGdxi,[4,1,2,3])),2)),eye(length(bhat_si))),1),2)); % 1/2*Tr(invG*dG)
            end
            
            if nderiv >= 3
                % second order derivatives
                
                ddphidbdb = Model.exp{s}.ddphidbdb(beta,bhat_si);
                ddphidbdbeta = Model.exp{s}.ddphidbdbeta(beta,bhat_si);
                
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
                
                if(Model.integration)
                    % laplace approximation
                    invG = pinv(G);
                    
                    ddGdbetadbeta = pdpdGpdbetapdbeta + 2*chainrule(permute(pddGdbpdbeta,[1,2,4,3]),dbhat_sidbeta) ...
                        + chainrule_ddxdydy_dydz(ddGdbdb,dbhat_sidbeta) + chainrule(dGdb,ddbhat_sidbetadbeta);
                    
                    ddGddeltaddelta = pdpdGpddeltapddelta + 2*chainrule(permute(pddGdbpddelta,[1,2,4,3]),dbhat_siddelta) ...
                        + chainrule_ddxdydy_dydz(ddGdbdb,dbhat_siddelta) + chainrule(dGdb,ddbhat_siddeltaddelta);
                    
                    ddGdbetaddelta = pdpdGpdbetapddelta + chainrule(permute(pddGdbpdbeta,[1,2,4,3]),dbhat_siddelta) ...
                        + permute(chainrule(permute(pddGdbpddelta,[1,2,4,3]),dbhat_sidbeta),[1,2,4,3]) ...
                        + chainrule_ddxdydy_dydz_dydv(ddGdbdb,dbhat_sidbeta,dbhat_siddelta) ...
                        + chainrule(dGdb,ddbhat_sidbetaddelta);
                    
                    ddGdxidxi = chainrule_ddxdydy_dydz(ddGdbetadbeta,dbetadxi) + chainrule(dGdbeta,ddbetadxidxi) ...
                        + chainrule_ddxdydy_dydz(ddGddeltaddelta,ddeltadxi) + chainrule(dGddelta,dddeltadxidxi) ...
                        + 2*chainrule_ddxdydy_dydz_dydv(ddGdbetaddelta,dbetadxi,ddeltadxi);
                    
                    dinvGdxi = squeeze(sum(bsxfun(@times,invG,permute(squeeze(sum(bsxfun(@times,permute(dGdxi,[4,1,2,3]),invG),2)),[4,1,2,3])),2));
                    
                    ddlogLi_Idxidxi(:,:,i) = - 0.5*squeeze(sum(sum(squeeze(bsxfun(@times,sum(bsxfun(@times,permute(dinvGdxi,[1,2,4,3]),permute(dGdxi,[4,1,2,5,3])),2)+sum(bsxfun(@times,invG,permute(ddGdxidxi,[5,1,2,3,4])),2),permute(eye(length(bhat_si)),[1,3,2]))),1),2)); % 1/2*Tr(dinvG*dg + invG*ddG)
                end
                
            end
            
        end
        
        Y_si_tmp{i} = Y_si;
        T_si_tmp{i} = T_si;
        R_si_tmp{i} = R_si;
        
        % Visulization
        if options.plot
            
            % Visualisation of single cell parameters
            if(isempty(fp))
                if(isfield(Model.exp{s},'title'))
                    if(ischar(Model.exp{s}.title))
                        fp(s) = figure('Name',Model.exp{s}.title);
                    else
                        fp(s) = figure;
                    end
                else
                    fp(s) = figure;
                end
            else
                if(length(fp)<s)
                    if(isfield(Model.exp{s},'title'))
                        if(ischar(Model.exp{s}.title))
                            fp(s) = figure('Name',Model.exp{s}.title);
                        else
                            fp(s) = figure;
                        end
                    else
                        fp(s) = figure;
                    end
                elseif(isempty(fp(s)))
                    if(isfield(Model.exp{s},'title'))
                        if(ischar(Model.exp{s}.title))
                            fp(s) = figure('Name',Model.exp{s}.title);
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
                %xlabel(char(Model.sym.b(Model.exp{s}.ind_b(j))));
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
                if(isfield(Model.exp{s},'title'))
                    if(ischar(Model.exp{s}.title))
                        fl(s) = figure('Name',Model.exp{s}.title);
                    else
                        fl(s) = figure;
                    end
                else
                    fl(s) = figure;
                end
            else
                if(length(fl)<s)
                    if(isfield(Model.exp{s},'title'))
                        if(ischar(Model.exp{s}.title))
                            fl(s) = figure('Name',Model.exp{s}.title);
                        else
                            fl(s) = figure;
                        end
                    else
                        fl(s) = figure;
                    end
                elseif(isempty(fl(s)))
                    if(isfield(Model.exp{s},'title'))
                        if(ischar(Model.exp{s}.title))
                            fl(s) = figure('Name',Model.exp{s}.title);
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
            if(Model.integration)
                bar([logLi_D;logLi_T;logLi_b;logLi_I;repmat(logL_s,[1,size(Data{s}.SCTL.Y,3)])],'stacked')
                set(gca,'XTickLabel',{'Data','Event','Par','Int','Pen'})
            else
                bar([logLi_D;logLi_T;logLi_b;repmat(logL_s,[1,size(Data{s}.SCTL.Y,3)])],'stacked')
                set(gca,'XTickLabel',{'Data','Event','Par','Pen'})
            end
            ylabel('log-likelihood')
            title('likelihood contribution')
            
            
            
            % Visualisation of data and fit
            Model.exp{s}.plot(Data{s},Sim_SCTL,s);
        end
        
    end
end
