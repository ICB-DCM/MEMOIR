% function [logL,dlogLdxi,ddlogLdxi2] = logL_CE_w_grad_2(xi,Data,Model,options,extract_flag)
function varargout = logL_CE_w_grad_2(varargin)
    
    %% Load
    persistent tau
    persistent P_old
    persistent logL_old
    
    if isempty(tau)
        tau = clock;
    end
    if isempty(logL_old)
        logL_old = -inf;
    end
    
    %% Initialization
    xi = varargin{1};
    Data = varargin{2};
    Model = varargin{3};
    
    % Options
    options.tau_update = 0;
    options.plot = 1;
    if nargin == 4
        options = setdefault(varargin{4},options);
    end
    
    if nargin >= 5
        extract_flag = varargin{5};
    else
        extract_flag = false;
    end
    
    % Plot options
    if (etime(clock,tau) > options.tau_update) && (options.plot == 1)
        options.plot = 30;
        tau = clock;
    else
        options.plot = 0;
    end
    
    %% Evaluation of likelihood function
    % Initialization
    logL = 0;
    if nargout >= 2
        dlogLdxi = zeros(length(xi),1);
        if nargout >= 3
            ddlogLdxidxi = zeros(length(xi));
        end
    end
    
    % Data types
    % - Single Cell Time-Lapse (SCTL)
    % - Population Average (PA)
    % - Single cell SnapsHot (SCSH)
    data_type = {'SCTL','PA','SCSH'};
    
    % Loop: Experimental conditions
    for s = 1:length(Data)
        
        %% Assignment of global variables
        type_D = Model.type_D;
        
        %% Construct fixed effects and covariance matrix
        %         xi_e = xi;
        %         xi_e(6) = xi(6) + 1e-3;
        beta = Model.exp{s}.beta(xi);
        %         beta_e = Model.exp{s}.beta(xi_e);
        delta = Model.exp{s}.delta(xi);
        %         delta_e = Model.exp{s}.delta(xi_e);
        
        [D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);
        %         [D_e,invD_e,dDddelta_e,dinvDddelta_e,ddDddeltaddelta_e,ddinvDddeltaddelta_e] = xi2D(delta_e,type_D);
        
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,1,3)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,3,5)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,2,4)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,4,6)
        
        
        %% Construction of time vector
        t_s = [];
        for dtype = 1:length(data_type)
            if isfield(Data{s},data_type{dtype})
                t_s = union(eval(['Data{s}.' data_type{dtype} '.time']),t_s);
            end
        end
        
        %% Single cell time-lapse data - Individuals
        if isfield(Data{s},'SCTL')
            % Reset values of Data likelihood and parameter likelihood
            logL_D = 0;
            logL_b = 0;
            logL_I = 0;
            
            % Evaluation of time index set
            [~,ind_t] = ismember(Data{s}.SCTL.time,t_s);
            
            % Initialization
            Sim_SCTL = nan(size(Data{s}.SCTL.Y));
            
            % Loop: Indiviudal cells
            for i = 1:size(Data{s}.SCTL.Y,3)
                % Load single-cell data
                Ym_si = Data{s}.SCTL.Y(ind_t,:,i);
                
                ind = find(~isnan(Ym_si));
                
                if logL_old == -inf
                    bhat_si0 = zeros(length(Model.exp{s}.sym.b),1);
                else
                    bhat_si0 = P_old{s}.SCTL.bhat(:,i);
                end
                
                % Higher order derivatives of the objective function for single cell parameters
                % here G is the Hessian of the objective function and bhat_si is the optimum of the objective function
                % with respect to b
                % F_diff and b_diff determine how many derivatives of the objective function and of the optimum need to
                % be computed
                switch(nargout)
                    case 0
                        F_diff = 0;
                        b_diff = 0;
                        [bhat_si,G] ...
                            = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s);
                    case 1
                        F_diff = 0;
                        b_diff = 0;
                        [bhat_si,G] ...
                            = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s);
                    case 2
                        b_diff = 1;
                        if(Model.integration)
                            F_diff = 3;
                            
                            [bhat_si,dbhat_sidbeta,dbhat_siddelta,...
                                G,dGdb,pdGpdbeta,pdGpddelta] ...
                                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s);
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,D,dDddelta,ddDddeltaddelta,invD,dinvDddelta,ddinvDddeltaddelta,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,1,2)
                        else
                            F_diff = 2;
                            [bhat_si,dbhat_sidbeta,dbhat_siddelta,...
                                G] ...
                                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s);
                            %                             [bhat_si_e,dbhat_sidbeta_e,dbhat_siddelta_e,G_e] ...
                            %                                 = optimize_SCTL_si(Model,Data,bhat_si0,beta_e,D_e,dDddelta_e,ddDddeltaddelta_e,invD_e,dinvDddelta_e,ddinvDddeltaddelta_e,t_s,Ym_si,ind,F_diff,b_diff,s);
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,1,2)
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,1,3)
                        end
                    case 3
                        if(Model.integration)
                            F_diff = 4;
                            b_diff = 2;
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,D,dDddelta,ddDddeltaddelta,invD,dinvDddelta,ddinvDddeltaddelta,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,2,4)
                            [bhat_si,dbhat_sidbeta,dbhat_siddelta,ddbhat_sidbetadbeta,ddbhat_sidbetaddelta,ddbhat_siddeltaddelta,...
                                G,dGdb,pdGpdbeta,pdGpddelta,ddGdbdb,pddGdbpdbeta,pdpdGpdbetapdbeta,pddGdbpddelta,pdpdGpddeltapddelta,...
                                pdpdGpdbetapddelta] ...
                                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s);
                        else
                            F_diff = 3;
                            b_diff = 2;
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,D,dDddelta,ddDddeltaddelta,invD,dinvDddelta,ddinvDddeltaddelta,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,2,4)
                            [bhat_si,dbhat_sidbeta,dbhat_siddelta,ddbhat_sidbetadbeta,ddbhat_sidbetaddelta,ddbhat_siddeltaddelta,...
                                G,dGdb,pdGpdbeta] ...
                                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,ind,F_diff,b_diff,s);
                        end
                end
                
                % Store bhat
                P{s}.SCTL.bhat(:,i) = bhat_si;
                
                % Construct single-cell parameter
                phi_si = Model.exp{s}.phi(beta,bhat_si);
                
                % Simulate model
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(phi_si,@(phi) Model.exp{s}.model(t_s,phi,Data{s}.condition),1e-4,4,6)
                if nargout <2
                    [~,~,~,Y] = Model.exp{s}.model(t_s,phi_si,Data{s}.condition);
                elseif nargout <3
                    [~,~,~,Y,~,dYdphi] = Model.exp{s}.model(t_s,phi_si,Data{s}.condition);
                else
                    if(Model.integration)
                        [~,~,~,Y,~,dYdphi,~,ddYdphidphi] = Model.exp{s}.model(t_s,phi_si,Data{s}.condition);
                    else
                        [~,~,~,Y,~,dYdphi] = Model.exp{s}.model(t_s,phi_si,Data{s}.condition);
                    end
                end
                Y_si = Y(ind_t,:);
                
                % Apply indexing to derivatives of Y
                if nargout >= 2
                    dY_sidphi = zeros(length(ind),length(phi_si));
                    if nargout >= 3
                        ddY_sidphidphi = zeros(length(ind),length(phi_si),length(phi_si));
                    end
                    for k = 1:length(phi_si)
                        temp = dYdphi(:,:,k);
                        dY_sidphi(:,k) = temp(ind);
                        if Model.integration
                            for j = 1:length(phi_si)
                                temp = ddYdphidphi(:,:,j,k);
                                ddY_sidphidphi(:,j,k) = temp(ind);
                            end
                        end
                    end
                end
                
                % Construct sigma
                sigma = Model.exp{s}.sigma(phi_si);
                
                % Adapt sigma to size of data
                if(size(sigma,1) == size(Ym_si,1))
                    if(size(sigma,2) == 1)
                        Sigma_si = repmat(sigma,[1,size(Ym_si,2)]);
                    elseif(size(sigma,2) == size(Ym_si,2))
                        Sigma_si = sigma;
                    else
                        error('Incompatible size of sigma parametrisation!')
                    end
                elseif(size(sigma,2) == size(Ym_si,2))
                    if(size(sigma,1) == 1)
                        Sigma_si = repmat(sigma,[size(Ym_si,1),1]);
                    else
                        error('Incompatible size of sigma parametrisation!')
                    end
                elseif(and(size(sigma,1)==1,size(sigma,2)==1))
                    Sigma_si = repmat(sigma,size(Ym_si));
                else
                    error('Incompatible size of sigma parametrisation!')
                end
                
                % Evaluation of likelihood and likelihood gradient
                
                % this is part accounts for the noise model
                % J_D = log(p(Y(b,beta)|D))
                switch(Model.exp{s}.noise_model)
                    case 'normal'
                        switch(nargout)
                            case 0
                                J_D = normal_noise(Y_si,Ym_si,Sigma_si,ind);
                            case 1
                                J_D = normal_noise(Y_si,Ym_si,Sigma_si,ind);
                            case 2
                                [J_D,dJ_DdY,dJ_DdSigma] = normal_noise(Y_si,Ym_si,Sigma_si,ind);
                                %                                 ee = zeros(size(Y_si));
                                %                                 ee(16) = 1e-10;
                                %                                 [J_D_ee] = normal_noise(Y_si,Ym_si,Sigma_si_e,ind);
                                %                                 [J_D_e,dJ_DdY_e,dJ_DdSigma_e] = normal_noise(Y_si_e,Ym_si,Sigma_si_e,ind);
                            case 3
                                [J_D,dJ_DdY,dJ_DdSigma,ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = normal_noise(Y_si,Ym_si,Sigma_si,ind);
                        end
                        
                    case 'lognormal'
                        switch(nargout)
                            case 0
                                J_D = lognormal_noise(Y_si,Ym_si,Sigma_si,ind);
                            case 1
                                J_D = lognormal_noise(Y_si,Ym_si,Sigma_si,ind);
                            case 2
                                [J_D,dJ_DdY,dJ_DdSigma] = lognormal_noise(Y_si,Ym_si,Sigma_si,ind);
                            case 3
                                [J_D,dJ_DdY,dJ_DdSigma,ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = lognormal_noise(Y_si,Ym_si,Sigma_si,ind);
                        end
                end
                
                % this part accounts for the parameter model
                % J_b = log(p(b_si|delta))
                switch(Model.exp{s}.parameter_model)
                    case 'normal'
                        switch(nargout)
                            case 0
                                J_b = normal_param(bhat_si,delta,type_D);
                            case 1
                                J_b = normal_param(bhat_si,delta,type_D);
                            case 2
                                [J_b,dJ_bdb,pdJ_bpddelta]= normal_param(bhat_si,delta,type_D);
                                %                                 [J_b_e,~,~]= normal_param(bhat_si_e,D_e,dDddelta_e,ddDddeltaddelta_e,invD_e,dinvDddelta_e,ddinvDddeltaddelta_e);
                            case 3
                                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(phi_si,@(b) normal_param(b,delta,type_D),1e-4,1,2)
                                [J_b,dJ_bdb,pdJ_bpddelta,ddJ_bdbdb,dpdJ_bdbpddelta,pdpdJ_bpddeltapddelta]= normal_param(bhat_si,delta,type_D);
                        end
                    case 'lognormal'
                        switch(nargout)
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
                
                % Adapt dsigmadphi to size of data
                if nargout >= 2
                    dsigmadphi = Model.exp{s}.dsigmadphi(phi_si);
                    dSigmadphi = zeros(length(ind),length(phi_si));
                    
                    if(size(dsigmadphi,1) == size(Ym_si,1))
                        if(size(dsigmadphi,2) == 1)
                            dSdphi = repmat(dsigmadphi,[1,size(Ym_si,2),1]);
                        elseif(size(dsigmadphi,2) == size(Ym_si,2))
                            dSdphi = dsigmadphi;
                        end
                    elseif(size(dsigmadphi,2) == size(Ym_si,2))
                        if(size(dsigmadphi,1) == 1)
                            dSdphi = repmat(dsigmadphi,[size(Ym_si,1),1,1]);
                        end
                    elseif(and(size(dsigmadphi,1)==1,size(dsigmadphi,2)==1))
                        dSdphi = repmat(dsigmadphi,[size(Ym_si),1]);
                    end
                    for k = 1:length(phi_si)
                        temp = dSdphi(:,:,k);
                        dSigmadphi(:,k) = temp(ind);
                    end
                end
                
                logL_D = logL_D - J_D;
                logL_b = logL_b - J_b;
                
                
                if(Model.integration)
                    % laplace approximation
                    logL_I = logL_I - 0.5*log(det(G));
                end
                
                if nargout >= 2
                    % first order derivatives
                    dbetadxi = Model.exp{s}.dbetadxi(xi);
                    ddeltadxi = Model.exp{s}.ddeltadxi(xi);
                    dphidb = Model.exp{s}.dphidb(beta,bhat_si);
                    pdphipdbeta  = Model.exp{s}.dphidbeta(beta,bhat_si);
                    
                    dphidbeta = chainrule_dxdy_dydz(dphidb,dbhat_sidbeta) + pdphipdbeta;
                    dphiddelta = chainrule_dxdy_dydz(dphidb,dbhat_siddelta);
                    dphidxi = chainrule_dxdy_dydz(dphidbeta,dbetadxi) + chainrule_dxdy_dydz(dphiddelta,ddeltadxi);
                    
                    dJ_Ddphi = chainrule_dxdy_dydz(dJ_DdY,dY_sidphi) + chainrule_dxdy_dydz(dJ_DdSigma,dSigmadphi) ;
                    dJ_Ddxi = chainrule_dxdy_dydz(dJ_Ddphi,dphidxi);
                    
                    dbdxi = chainrule_dxdy_dydz(dbhat_sidbeta,dbetadxi) + chainrule_dxdy_dydz(dbhat_siddelta,ddeltadxi);
                    P{s}.SCTL.dbdxi(:,:,i) = dbdxi;
                    
                    dJ_bdxi = chainrule_dxdy_dydz(dJ_bdb,dbdxi) + chainrule_dxdy_dydz(pdJ_bpddelta,ddeltadxi);
                    
                    dlogLdxi = dlogLdxi - transpose(dJ_Ddxi) - transpose(dJ_bdxi);
                    
                    if(Model.integration)
                        % laplace approximation
                        invG = pinv(G);
                        
                        % take care when nelem(b) == 1 ... (ones(1,1,1) ==
                        % ones(1,1) so dGdb will be missing one dimension!)
                        if(numel(bhat_si)==1)
                            dGdbeta = pdGpdbeta + permute(dGdb*dbhat_sidbeta,[3,1,2]);
                            dGddelta = pdGpddelta + permute(dGdb*dbhat_siddelta,[3,1,2]);
                            dGdxi = chainrule_dxdy_dydz(dGdbeta,dbetadxi) + chainrule_dxdy_dydz(dGddelta,ddeltadxi);
                        else
                            dGdbeta = pdGpdbeta + chainrule_dxdy_dydz(dGdb,dbhat_sidbeta);
                            dGddelta = pdGpddelta + chainrule_dxdy_dydz(dGdb,dbhat_siddelta);
                            dGdxi = chainrule_dxdy_dydz(dGdbeta,dbetadxi) + chainrule_dxdy_dydz(dGddelta,ddeltadxi);
                        end
                        
                        
                        
                        dlogLdxi = dlogLdxi ...
                            - 0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,invG,permute(dGdxi,[4,1,2,3])),2)),eye(length(bhat_si))),1),2)); % 1/2*Tr(invG*dG)
                    end
                    
                    if nargout >= 3
                        % second order derivatives
                        ddbetadxidxi = Model.exp{s}.ddbetadxidxi(xi);
                        dddeltadxidxi = Model.exp{s}.dddeltadxidxi(xi);
                        
                        ddphidbdb = Model.exp{s}.ddphidbdb(beta,bhat_si);
                        ddphidbdbeta = Model.exp{s}.ddphidbdbeta(beta,bhat_si);
                        
                        ddphidbetadbeta = chainrule_dxdy_dydz(dphidb,ddbhat_sidbetadbeta) + chainrule_ddxdydy_dydz(ddphidbdb,dbhat_sidbeta) ...
                            + permute(chainrule_dxdy_dydz(permute(ddphidbdbeta,[1,3,2]),dbhat_sidbeta),[1,3,2]);
                        ddphidbetaddelta = chainrule_dxdy_dydz(dphidb,ddbhat_sidbetaddelta) + chainrule_ddxdydy_dydz_dydv(ddphidbdb,dbhat_sidbeta,dbhat_siddelta);
                        ddphiddeltaddelta = chainrule_dxdy_dydz(dphidb,ddbhat_siddeltaddelta) + chainrule_ddxdydy_dydz(ddphidbdb,dbhat_siddelta);
                        
                        ddphidxidxi = chainrule_dxdy_dydz(dphidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddphidbetadbeta,dbetadxi) ...
                            + chainrule_dxdy_dydz(dphiddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddphiddeltaddelta,ddeltadxi) ...
                            + chainrule_ddxdydy_dydz_dydv(ddphidbetaddelta,dbetadxi,ddeltadxi) + permute(chainrule_ddxdydy_dydz_dydv(ddphidbetaddelta,dbetadxi,ddeltadxi),[1,3,2]);
                        
                        ddJ_DdphidY = bsxfun(@times,ddJ_DdYdY,permute(dY_sidphi,[3,1,2])) + ...
                            bsxfun(@times,ddJ_DdYdSigma,permute(dSigmadphi,[3,1,2]));
                        
                        ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dY_sidphi,[3,1,2])) + ...
                            bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigmadphi,[3,1,2]));
                        
                        ddJ_Ddphidphi = chainrule_dxdy_dydz(dJ_DdY,ddY_sidphidphi) ...
                            + squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(dY_sidphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddJ_DdphidSigma,permute(dSigmadphi,[3,1,4,2])),2));
                        
                        ddJ_Ddxidxi = chainrule_dxdy_dydz(dJ_Ddphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Ddphidphi,dphidxi);
                        
                        ddbdxidxi = chainrule_dxdy_dydz(dbhat_sidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddbhat_sidbetadbeta,dbetadxi) ...
                            + chainrule_dxdy_dydz(dbhat_siddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddbhat_siddeltaddelta,ddeltadxi) ...
                            + chainrule_ddxdydy_dydz_dydv(ddbhat_sidbetaddelta,dbetadxi,ddeltadxi) ...
                            + chainrule_ddxdydy_dydz_dydv(permute(ddbhat_sidbetaddelta,[1,3,2]),ddeltadxi,dbetadxi);
                        
                        P{s}.SCTL.ddbdxidxi(:,:,:,i) = ddbdxidxi;
                        
                        ddJ_bdxidxi = chainrule_dxdy_dydz(dJ_bdb,ddbdxidxi) + chainrule_ddxdydy_dydz(ddJ_bdbdb,dbdxi) ...
                            + chainrule_dxdy_dydz(pdJ_bpddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(pdpdJ_bpddeltapddelta,ddeltadxi) ...
                            + chainrule_ddxdydy_dydz_dydv(dpdJ_bdbpddelta,dbdxi,ddeltadxi) ...
                            + chainrule_ddxdydy_dydz_dydv(permute(dpdJ_bdbpddelta,[2,1]),ddeltadxi,dbdxi);
                        
                        ddlogLdxidxi = ddlogLdxidxi - ddJ_Ddxidxi - ddJ_bdxidxi;
                        
                        if(Model.integration)
                            % laplace approximation
                            invG = pinv(G);
                            
                            ddGdbetadbeta = pdpdGpdbetapdbeta + 2*chainrule_dxdy_dydz(permute(pddGdbpdbeta,[1,2,4,3]),dbhat_sidbeta) ...
                                + chainrule_ddxdydy_dydz(ddGdbdb,dbhat_sidbeta) + chainrule_dxdy_dydz(dGdb,ddbhat_sidbetadbeta);
                            
                            ddGddeltaddelta = pdpdGpddeltapddelta + 2*chainrule_dxdy_dydz(permute(pddGdbpddelta,[1,2,4,3]),dbhat_siddelta) ...
                                + chainrule_ddxdydy_dydz(ddGdbdb,dbhat_siddelta) + chainrule_dxdy_dydz(dGdb,ddbhat_siddeltaddelta);
                            
                            ddGdbetaddelta = pdpdGpdbetapddelta + chainrule_dxdy_dydz(permute(pddGdbpdbeta,[1,2,4,3]),dbhat_siddelta) ...
                                + permute(chainrule_dxdy_dydz(permute(pddGdbpddelta,[1,2,4,3]),dbhat_sidbeta),[1,2,4,3]) ...
                                + chainrule_ddxdydy_dydz_dydv(ddGdbdb,dbhat_sidbeta,dbhat_siddelta) ...
                                + chainrule_dxdy_dydz(dGdb,ddbhat_sidbetaddelta);
                            
                            ddGdxidxi = chainrule_ddxdydy_dydz(ddGdbetadbeta,dbetadxi) + chainrule_dxdy_dydz(dGdbeta,ddbetadxidxi) ...
                                + chainrule_ddxdydy_dydz(ddGddeltaddelta,ddeltadxi) + chainrule_dxdy_dydz(dGddelta,dddeltadxidxi) ...
                                + 2*chainrule_ddxdydy_dydz_dydv(ddGdbetaddelta,dbetadxi,ddeltadxi);
                            
                            dinvGdxi = squeeze(sum(bsxfun(@times,invG,permute(squeeze(sum(bsxfun(@times,permute(dGdxi,[4,1,2,3]),invG),2)),[4,1,2,3])),2));
                            
                            ddlogLdxidxi = ddlogLdxidxi ...
                                - 0.5*squeeze(sum(sum(squeeze(bsxfun(@times,sum(bsxfun(@times,permute(dinvGdxi,[1,2,4,3]),permute(dGdxi,[4,1,2,5,3])),2)+sum(bsxfun(@times,invG,permute(ddGdxidxi,[5,1,2,3,4])),2),permute(eye(length(bhat_si)),[1,3,2]))),1),2)); % 1/2*Tr(dinvG*dg + invG*ddG)
                        end
                        
                    end
                    
                end
                
                % Assignment
                Sim_SCTL(:,:,i) = Y(ind_t,:);
            end
            
            logL = logL + logL_D + logL_b + logL_I;
            
            if(Model.penalty)
                % parameter penalization terms
                % logL_s = log(p(mu_S,S_s|b_s,D))
                P{s}.SCTL.bhat(:,i) = bhat_si;
                if nargout<= 1
                    logL_s = penal_param(P{s}.SCTL.bhat,delta,type_D);
                elseif nargout<= 2
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(P{s}.SCTL.bhat,@(x) penal_param(x,delta,type_D),1e-4,1,2)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) penal_param(P{s}.SCTL.bhat,x,type_D),1e-4,1,3)
                    [logL_s,dlogL_sdb_s,dlogL_sddelta] = penal_param(P{s}.SCTL.bhat,delta,type_D);
                elseif nargout<= 3
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(P{s}.SCTL.bhat,@(x) penal_param(x,delta,type_D),1e-4,2,4)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) penal_param(P{s}.SCTL.bhat,x,type_D),1e-4,2,5)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) penal_param(P{s}.SCTL.bhat,x,type_D),1e-4,3,6)
                    [logL_s,dlogL_sdb_s,dlogL_sddelta,ddlogL_sdb_sdb_s,ddlogL_sdb_sddelta,ddlogL_sddeltaddelta] = penal_param(P{s}.SCTL.bhat,delta,type_D);
                end
                
                
                logL = logL + logL_s;
                
                if nargout >= 2
                    dlogL_sdxi = squeeze(sum(sum(bsxfun(@times,dlogL_sdb_s,permute(P{s}.SCTL.dbdxi,[1,3,2])),1),2)) + transpose(chainrule_dxdy_dydz(dlogL_sddelta,ddeltadxi));
                    dlogLdxi = dlogLdxi + dlogL_sdxi;
                    if nargout >= 3
                        ddlogL_sdxidxi = squeeze(sum(sum(bsxfun(@times,bsxfun(@times,ddlogL_sdb_sdb_s,permute(P{s}.SCTL.dbdxi,[1,3,2])),permute(P{s}.SCTL.dbdxi,[1,3,4,2])),1),2)) ...
                            + squeeze(sum(sum(bsxfun(@times,dlogL_sdb_s,permute(P{s}.SCTL.ddbdxidxi,[1,4,2,3])),1),2)) ...
                            + squeeze(sum(sum(sum(bsxfun(@times,bsxfun(@times,ddlogL_sdb_sddelta,permute(P{s}.SCTL.dbdxi,[1,3,4,2])),permute(ddeltadxi,[3,4,1,5,2])),1),2),3)) ...
                            + transpose(squeeze(sum(sum(sum(bsxfun(@times,bsxfun(@times,ddlogL_sdb_sddelta,permute(P{s}.SCTL.dbdxi,[1,3,4,2])),permute(ddeltadxi,[3,4,1,5,2])),1),2),3))) ...
                            + chainrule_ddxdydy_dydz(ddlogL_sddeltaddelta,ddeltadxi) ...
                            + chainrule_dxdy_dydz(dlogL_sddelta,dddeltadxidxi);
                        ddlogLdxidxi = ddlogLdxidxi + ddlogL_sdxidxi;
                    end
                    
                end
            end
            
            
            
            
            % Visulization
            if options.plot
                
                % Visualisation of single cell parameters
                figure(Model.exp{s}.fp)
                clf
                b_s = P{s}.SCTL.bhat;
                d_s = size(b_s,2);
                n_b = size(b_s,1);
                
                mu_s = 1/d_s*(sum(b_s,2));
                S_s = 1/d_s*squeeze(sum(bsxfun(@times,permute(b_s,[2,1]),permute(b_s,[2,3,1])),1)) + 1e-10*eye(n_b); % regularization
                for j = 1:n_b
                    subplot(ceil((n_b+1)/4),4,j+1)
                    xx = linspace(-5*sqrt(D(j,j)),5*sqrt(D(j,j)),100);
                    %nhist(P{s}.SCTL.bhat(j,:),'pdf','noerror');
                    hold on
                    plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','LineWidth',2)
                    plot(xx,normcdf(xx,0,sqrt(S_s(j,j))),'--r','LineWidth',2)
                    
                    
                    if(j==1)
                        
                    end
                    xlim([-5*sqrt(D(j,j)),5*sqrt(D(j,j))])
                    ylim([0,1.1])
                    box on
                end
                subplot(ceil(n_b+1/4),4,1,'Visible','off')
                hold on
                plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','Visible','off')
                plot(xx,normcdf(xx,0,sqrt(S_s(j,j))),'--r','Visible','off')
                
                
                legend('cdf of single cell Parameters','cdf of populaton Parameters')
                
                % Visualisation of likelihood contribution
                
                figure(Model.exp{s}.fl)
                if(Model.penalty)
                    if(Model.integration)
                        bar([logL_D,logL_b,logL_I,logL_s])
                        set(gca,'XTickLabel',{'Data','Par','Int','Pen'})
                    else
                        bar([logL_D,logL_b,logL_s])
                        set(gca,'XTickLabel',{'Data','Par','Pen'})
                    end
                else
                    if(Model.integration)
                        bar([logL_D,logL_b,logL_I])
                        set(gca,'XTickLabel',{'Data','Par','Int'})
                    else
                        bar([logL_D,logL_b])
                        set(gca,'XTickLabel',{'Data','Par'})
                    end
                end
                ylabel('log-likelihood')
                title('likelihood contribution')
                
                
                
                % Visualisation of data and fit
                Model.exp{s}.plot(Data{s},Sim_SCTL,Model.exp{s}.fh);
            end
        end
        
        %% Single cell time-lapse data - Statistics
        if isfield(Data{s},'SCTLstat')
            % Simulation using sigma points
            if nargout == 1
                [~,~,~,mz_SP,Cz_SP,B_SP,~] = ...
                    getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition),beta,D,Model.exp{s});
            else
                dbetadxi = Model.exp{s}.dbetadxi(xi);
                ddeltadxi = Model.exp{s}.ddeltadxi(xi);
                [~,~,~,mz_SP,Cz_SP,B_SP,~,~,~,~,dmz_SPdxi,dCz_SPdxi,~,~] = ...
                    getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition),beta,D,Model.exp{s},dDddelta,ddeltadxi,dbetadxi);
            end
            
            % Evaluation of likelihood, likelihood gradient and hessian
            
            % Mean
            logL_mz = - 0.5*sum(sum(((Data{s}.SCTLstat.mz - mz_SP)./Data{s}.SCTLstat.Sigma_mz).^2,1),2);
            if nargout >= 2
                dlogL_mzdxi = permute(sum(bsxfun(@times,(Data{s}.SCTLstat.mz - mz_SP)./Data{s}.SCTLstat.Sigma_mz.^2,dmz_SPdxi),1),[2,1]);
                if nargout >= 3
                    wdmz_SP = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_mz,dmz_SPdxi);
                    %                     wdmz_SP = reshape(wdmz_SP,[numel(mz_SP),size(dmz_SPdxi,3)]);
                    ddlogL_mzdxi2 = -wdmz_SP'*wdmz_SP;
                end
            end
            
            
            % Covariance
            logL_Cz = - 0.5*sum(sum(sum(((Data{s}.SCTLstat.Cz - Cz_SP)./Data{s}.SCTLstat.Sigma_Cz).^2,1),2),3);
            if nargout >= 2
                dlogL_Czdxi = squeeze(sum(sum(bsxfun(@times,(Data{s}.SCTLstat.Cz - Cz_SP)./Data{s}.SCTLstat.Sigma_Cz.^2,dCz_SPdxi),1),2));
                if nargout >= 3
                    wdCz_SP = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_Cz,dCz_SPdxi);
                    wdCz_SP = reshape(wdCz_SP,[numel(Cz_SP),size(dCz_SPdxi,3)]);
                    ddlogL_Czdxi2 = -wdCz_SP'*wdCz_SP;
                end
            end
            
            % Summation
            logL = logL + logL_mz + logL_Cz;
            if nargout >=2
                dlogLdxi = dlogLdxi + dlogL_mzdxi + dlogL_Czdxi;
                if nargout >=3
                    ddlogLdxidxi = ddlogLdxidxi + ddlogL_mzdxi2 + ddlogL_Czdxi2;
                end
            end
            
            % Visulization
            if options.plot
                Sim_SCTLstat.mz = mz_SP;
                Sim_SCTLstat.Cz = Cz_SP;
                Model.exp{s}.plot(Data{s},Sim_SCTLstat,Model.exp{s}.fh);
            end
            
        end
        
        %% Single cell snapshot data
        if isfield(Data{s},'SCSH')
            % Simulation using sigma points
            if nargout == 1
                [m_SP,C_SP,~,~,~,B_SP,~] = ...
                    getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCSH.time,phi,Data{s}.condition),beta,D,Model.exp{s});
            else
                dbetadxi = Model.exp{s}.dbetadxi(xi);
                ddeltadxi = Model.exp{s}.ddeltadxi(xi);
                [m_SP,C_SP,~,~,~,B_SP,~,dm_SP,dC_SP,~,~,~,~,~] = ...
                    getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCSH.time,phi,Data{s}.condition),beta,D,Model.exp{s},dDddelta,ddeltadxi,dbetadxi);
            end
            
            % Evaluation of likelihood, likelihood gradient and hessian
            % Mean
            logL_m = - 0.5*sum(sum(((Data{s}.SCSH.m - m_SP)./Data{s}.SCSH.Sigma_m).^2,1),2);
            if nargout >= 2
                dlogL_mdxi = squeeze(sum(sum(bsxfun(@times,(Data{s}.SCSH.m - m_SP)./Data{s}.SCSH.Sigma_m.^2,dm_SP),1),2));
                if narogut >= 3
                    wdm_SP = bsxfun(@times,1./Data{s}.SCSH.Sigma_m,dm_SP);
                    wdm_SP = reshape(wdm_SP,[numel(m_SP),size(dm_SP,3)]);
                    ddlogL_mdxi2 = -wdm_SP'*wdm_SP;
                end
            end
            
            % Covariance
            logL_C = - 0.5*sum(sum(sum(((Data{s}.SCSH.C - C_SP)./Data{s}.SCSH.Sigma_C).^2,1),2),3);
            if nargout >= 2
                dlogL_Cdxi = squeeze(sum(sum(sum(bsxfun(@times,(Data{s}.SCSH.C - C_SP)./Data{s}.SCSH.Sigma_C.^2,dC_SP),1),2),3));
                if nargout >= 3
                    wdC_SP = bsxfun(@times,1./Data{s}.SCSH.Sigma_C,dC_SP);
                    wdC_SP = reshape(wdC_SP,[numel(C_SP),size(dC_SP,4)]);
                    ddlogL_Cdxi2 = -wdC_SP'*wdC_SP;
                end
            end
            
            % Summation
            logL = logL + logL_m + logL_C;
            if nargout >= 2
                dlogLdxi = dlogLdxi + dlogL_mdxi + dlogL_Cdxi;
                if nargout >= 3
                    ddlogLdxidxi = ddlogLdxidxi + ddlogL_mdxi2 + ddlogL_Cdxi2;
                end
            end
            
            % Visulization
            if options.plot
                Sim_SCSH.m = m_SP;
                Sim_SCSH.C = C_SP;
                Model.exp{s}.plot(Data{s},Sim_SCSH,Model.exp{s}.fh);
            end
            
        end
        
        %% Population average data
        if isfield(Data{s},'PA')
            
            % Simulation using sigma points
            if nargout == 1
                [m_SP,~,~,~,~,B_SP,~] = ...
                    getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.PA.time,phi,Data{s}.condition),beta,D,Model.exp{s});
            else
                dbetadxi = Model.exp{s}.dbetadxi(xi);
                ddeltadxi = Model.exp{s}.ddeltadxi(xi);
                [m_SP,~,~,~,~,B_SP,~,dm_SP,~,~,~,~,~,~] = ...
                    getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.PA.time,phi,Data{s}.condition),beta,D,Model.exp{s},dDddelta,ddeltadxi,dbetadxi);
            end
            
            % Post-processing of population average data
            if isfield(Model.exp{s},'PA_post_processing')
                if(nargout==1) 
                    dm_SP = zeros([size(m_SP) size(xi,1)]);
                end
                [m_SP,dm_SP] = Model.exp{s}.PA_post_processing(m_SP,dm_SP);
            end
            
            
            % Evaluation of likelihood, likelihood gradient and hessian
            logL_m = - 0.5*sum(sum(((Data{s}.PA.m - m_SP)./Data{s}.PA.Sigma_m).^2,1),2);
            if nargout >= 2
                dlogL_mdxi = squeeze(sum(sum(bsxfun(@times,(Data{s}.PA.m - m_SP)./Data{s}.PA.Sigma_m.^2,dm_SP),1),2));
                if nargout >= 3
                    wdm_SP = bsxfun(@times,1./Data{s}.PA.Sigma_m,dm_SP);
                    wdm_SP = reshape(wdm_SP,[numel(m_SP),size(dm_SP,3)]);
                    ddlogL_mdxi2 = -wdm_SP'*wdm_SP;
                end
            end
            
            % Summation
            logL = logL + logL_m;
            if nargout >= 2
                dlogLdxi = dlogLdxi + dlogL_mdxi;
                if nargout >= 3
                    ddlogLdxidxi = ddlogLdxidxi + ddlogL_mdxi2;
                end
            end
            
            % Visulization
            if options.plot
                Sim_PA.m = m_SP;
                Model.exp{s}.plot(Data{s},Sim_PA,Model.exp{s}.fh);
            end
        end
        
    end
    
    %% Output
    
    if extract_flag
        if isfield(Data{s},'SCTL')
            varargout{1} = P;
        elseif any([isfield(Data{s},'SCTLstat'),isfield(Data{s},'PA'),isfield(Data{s},'SCSH')])
            varargout{1} = B_SP;
        end
        return
    end
    
    if nargout >= 1
        % One output
        varargout{1} =  logL;
        if nargout >= 2
            % Two outputs
            varargout{2} =  dlogLdxi;
            if nargout >= 3
                % Two outputs
                varargout{3} =  ddlogLdxidxi;
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTION FOR OPTIMIZATION OF SINGLE CELL PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [bhat_si,dbhat_sidbeta,dbhat_siddelta,ddbhat_sidbetadbeta,ddbhat_sidbetaddelta,ddbhat_siddeltaddelta,dGdbeta,dGddelta,ddGdbetadbeta,ddGdbetaddelta,ddGddeltaddelta] = optimize_SCTL_si(Model,Data,beta,invD,t_s,Ym_si,Sigma_si,ind,n_b,n_xi,F_diff,b_diff);
function varargout = optimize_SCTL_si(Model,Data,bhat_0,beta,delta,type_D,t,Ym,ind,F_diff,b_diff,s)
    options_fmincon = optimset('algorithm','trust-region-reflective',...
        'display','off',...
        'GradObj','on',...
        'MaxIter',1000,... % 1000%'display','iter',...
        'TolFun',1e-6,...
        'TolX',1e-6,...
        'Hessian','user-supplied');
    
    [bhat,~,~,~,~,~,~] = fmincon(...
        @(b) objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,ind,s),...
        bhat_0,[],[],[],[],-10*ones(length(bhat_0),1),10*ones(length(bhat_0),1),[],options_fmincon);
    
    
    switch(b_diff)
        case 0
            % Higher order derivatives objective function
            [~,~,G] = objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s);
            varargout{1} = bhat;
            varargout{2} = G;
        case 1
            switch(F_diff)
                case 3
                    % Higher order derivatives objective function
                    [~,~,G,...
                        ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                        dGdb,pdGpdbeta,pdGpddelta] = ...
                        objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s);
                    % Gradient of optimal point
                    [bhat,dbhatdbeta,dbhatddelta] = bhat_SCTL_si(bhat,G,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta);
                    
                    varargout{1} = bhat;
                    
                    varargout{2} = dbhatdbeta;
                    varargout{3} = dbhatddelta;
                    
                    varargout{4} = G;
                    
                    varargout{5} = dGdb;
                    varargout{6} = pdGpdbeta;
                    varargout{7} = pdGpddelta;
                case 2
                    % Higher order derivatives objective function
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,ind,s),1e-4,2,3)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(beta,@(beta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s),1e-4,2,4)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(delta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s),1e-4,2,5)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,ind,s),1e-4,1,2)
                    [~,~,G,...
                        ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta] = ...
                        objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s);
                    % Gradient of optimal point
                    [bhat,dbhatdbeta,dbhatddelta] = bhat_SCTL_si(bhat,G,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta);
                    varargout{1} = bhat;
                    varargout{2} = dbhatdbeta;
                    varargout{3} = dbhatddelta;
                    varargout{4} = G;
            end
        case 2
            switch(F_diff)
                case 4
                    % Higher order derivatives objective function
                    [~,~,G,...
                        ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                        dGdb,pdGpdbeta,pdGpddelta,...
                        dddJdbdbetadbeta,dddJdbddeltaddelta,...
                        ddGdbdb,pddGdbpdbeta,pdpdGpdbetapdbeta,pddGdbpddelta,pdpdGpddeltapddelta,pdpdGpdbetapddelta] = ...
                        objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s);
                    % Gradient + Hessian of optimal point
                    [bhat,dbhatdbeta,dbhatddelta,ddbhatdbetadbeta,ddbhatdbetaddelta,ddbhatddeltaddelta] = bhat_SCTL_si(bhat,G,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,dGdb,pdGpdbeta,pdGpddelta,dddJdbdbetadbeta,dddJdbddeltaddelta);
                    
                    
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
                    [~,~,G,...
                        ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                        dGdb,pdGpdbeta,pdGpddelta,...
                        dddJdbdbetadbeta,dddJdbddeltaddelta] = ...
                        objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,ind,s);
                    % Gradient + Hessian of optimal point
                    [bhat,dbhatdbeta,dbhatddelta,ddbhatdbetadbeta,ddbhatdbetaddelta,ddbhatddeltaddelta] = bhat_SCTL_si(bhat,G,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,dGdb,pdGpdbeta,pdGpddelta,dddJdbdbetadbeta,dddJdbddeltaddelta);
                    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTION FOR DERIVATIVE OF OPTIMA W.R.T HYPERPARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [bhat,dbhatdbeta,dbhatddelta,ddbhatdbetadbeta,ddbhatdbetaddelta,ddbhatddeltaddelta] = bhat_SCTL_si(bhat,G,pdGpdbeta,pdGpddelta);
function varargout = bhat_SCTL_si(bhat,G,ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,dGdb,pdGpdbeta,pdGpddelta,dddJdbdbetadbeta,dddJdbddeltaddelta)
    
    invG = pinv(G);
    
    if nargout>=2
        
        dbhatdbeta = -invG*squeeze(ddJdbdbeta);
        dbhatddelta = -invG*squeeze(ddJdbddelta);
        
        if nargout >= 4
            
            invGdGdbeta = chainrule_dxdy_dydz(invG,pdGpdbeta);
            invGdGddelta = chainrule_dxdy_dydz(invG,pdGpddelta);
            
            ddbhatdbetadbeta  = permute(sum(bsxfun(@times,permute(invGdGdbeta,[1,2,4,3]),permute(dbhatdbeta,[3,1,2])),2),[1,3,4,2]) - chainrule_dxdy_dydz(invG,dddJdbdbetadbeta);
            
            ddbhatddeltaddelta = permute(sum(bsxfun(@times,permute(invGdGddelta,[1,2,4,3]),permute(dbhatddelta,[3,1,2])),2),[1,3,4,2]) - chainrule_dxdy_dydz(invG,dddJdbddeltaddelta);
            
            ddbhatdbetaddelta = permute(sum(bsxfun(@times,permute(invGdGddelta,[1,2,4,3]),permute(dbhatdbeta,[3,1,2])),2),[1,3,4,2]);%- chainrule_dxdy_dydz(invG,dddJdbdbetaddelta) = 0
            
        end
    end
    
    
    
    if nargout <= 1
        varargout{1} =  bhat;
    elseif nargout <= 3
        varargout{1} =  bhat;
        varargout{2} =  dbhatdbeta;
        varargout{3} =  dbhatddelta;
    elseif nargout <= 6
        varargout{1} =  bhat;
        varargout{2} =  dbhatdbeta;
        varargout{3} =  dbhatddelta;
        varargout{4} =  ddbhatdbetadbeta;
        varargout{5} =  ddbhatdbetaddelta;
        varargout{6} =  ddbhatddeltaddelta;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OBJECTIVE FUNCTION FOR OPTIMIZATION OF SINGLE CELL PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [~,~,G,ddJdbdbeta,ddJdbdelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,Gdb,dGdsigma,ddGdbdb,ddGdbdsigma,ddGdsigmadsigma] = objective_SCTL_s1_full(Model,beta,bhat,Data{s}.condition,invD,dinvDddelta,ddinvDddeltaddelta,t_s,Ym,Sigma,ind,s);
% function varargout = objective_SCTL_s1(Model,beta,b,kappa,D,dDddelta,ddDddeltaddelta,invD,dinvDddelta,ddinvDddeltaddelta,t,Ym,ind,s)

function varargout = objective_SCTL_s1(Model,beta,b,kappa,delta,type_D,t,Ym,ind,s)

    [D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);
    
    % Single-cell parameters
    phi = Model.exp{s}.phi(beta,b);
    
    sigma = Model.exp{s}.sigma(phi);
    
    if(size(sigma,1) == size(Ym,1))
        if(size(sigma,2) == 1)
            Sigma = repmat(sigma,[1,size(Ym,2)]);
        elseif(size(sigma,2) == size(Ym,2))
            Sigma = sigma;
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == size(Ym,2))
        if(size(sigma,1) == 1)
            Sigma = repmat(sigma,[size(Ym,1),1]);
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        Sigma = repmat(sigma,size(Ym));
    else
        error('Incompatible size of sigma parametrisation!')
    end
    
    % Simulate model
    if nargout >= 3
        %[~,~,~,Y,~,dYdphi,~,ddYdphidphi] = Model.exp{s}.model(t,phi,kappa);
        [~,~,~,Y,~,dYdphi] = Model.exp{s}.model(t,phi,kappa);
        temp1 = dYdphi;
        %         temp2 = ddYdphidphi;
        dYdphi = zeros(length(ind),length(phi));
        ddYdphidphi = zeros(length(ind),length(phi),length(phi));
        for j = 1:length(phi)
            temp = temp1(:,:,j);
            dYdphi(:,j) = temp(ind);
            %             for k = 1:length(phi)
            %                 temp = temp2(:,:,j,k);
            %                 ddYdphidphi(:,j,k) = temp(ind);
            %             end
        end
    elseif nargout >=2
        [~,~,~,Y,~,dYdphi] = Model.exp{s}.model(t,phi,kappa);
        temp1 = dYdphi;
        dYdphi = zeros(length(ind),length(phi));
        for j = 1:length(phi)
            temp = temp1(:,:,j);
            dYdphi(:,j) = temp(ind);
        end
    else
        [~,~,~,Y] = Model.exp{s}.model(t,phi,kappa);
    end
    
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)Model.exp{s}.model(t,phi,kappa),1e-4,4,6)
    
    if nargout >= 2 % first order derivatives
        dsigmadphi = Model.exp{s}.dsigmadphi(phi);
        dSigmadphi = zeros(length(ind),length(phi));
        
        if(size(dsigmadphi,1) == size(Ym,1))
            if(size(dsigmadphi,2) == 1)
                dSdphi = repmat(dsigmadphi,[1,size(Ym,2),1]);
            elseif(size(dsigmadphi,2) == size(Ym,2))
                dSdphi = dsigmadphi;
            end
        elseif(size(dsigmadphi,2) == size(Ym,2))
            if(size(dsigmadphi,1) == 1)
                dSdphi = repmat(dsigmadphi,[size(Ym,1),1,1]);
            end
        elseif(and(size(dsigmadphi,1)==1,size(dsigmadphi,2)==1))
            dSdphi = repmat(dsigmadphi,[size(Ym),1]);
        end
        
        if nargout >= 3 % second order derivatives
            ddsigmadphidphi = Model.exp{s}.ddsigmadphidphi(phi);
            ddSigmadphidphi = zeros(length(ind),length(phi),length(phi));
            if(size(dsigmadphi,1) == size(Ym,1))
                if(size(dsigmadphi,2) == 1)
                    ddSdphidphi = repmat(ddsigmadphidphi,[1,size(Ym,2),1,1]);
                elseif(size(dsigmadphi,2) == size(Ym,2))
                    ddSdphidphi = ddsigmadphidphi ;
                end
            elseif(size(dsigmadphi,2) == size(Ym,2))
                if(size(dsigmadphi,1) == 1)
                    ddSdphidphi = repmat(ddsigmadphidphi,[size(Ym,1),1,1]);
                end
            elseif(and(size(dsigmadphi,1)==1,size(dsigmadphi,2)==1))
                ddSdphidphi = repmat(ddsigmadphidphi,[size(Ym),1,1]);
            end
        end
        
        if nargout >= 9 % third order derivatives
            dddsigmadphidphidphi = Model.exp{s}.dddsigmadphidphidphi(phi);
            dddSigmadphidphidphi = zeros(length(ind),length(phi),length(phi),length(phi));
            if(size(dsigmadphi,1) == size(Ym,1))
                if(size(dsigmadphi,2) == 1)
                    dddSdphidphidphi = repmat(dddsigmadphidphidphi,[1,size(Ym,2),1,1,1]);
                elseif(size(dsigmadphi,2) == size(Ym,2))
                    dddSdphidphidphi = dddsigmadphidphidphi;
                end
            elseif(size(dsigmadphi,2) == size(Ym,2))
                if(size(dsigmadphi,1) == 1)
                    dddSdphidphidphi = repmat(dddsigmadphidphidphi,[size(Ym,1),1,1,1]);
                end
            elseif(and(size(dsigmadphi,1)==1,size(dsigmadphi,2)==1))
                dddSdphidphidphi = repmat(dddsigmadphidphidphi,[size(Ym),1,1,1]);
            end
        end
        
        if nargout >= 14 % fourth order derivatives
            ddddsigmadphidphidphidphi = Model.exp{s}.ddddsigmadphidphidphidphi(phi);
            ddddSigmadphidphidphidphi = zeros(length(ind),length(phi),length(phi),length(phi),length(phi));
            if(size(dsigmadphi,1) == size(Ym,1))
                if(size(dsigmadphi,2) == 1)
                    ddddSdphidphidphidphi = repmat(ddddsigmadphidphidphidphi,[1,size(Ym,2),1,1,1]);
                elseif(size(dsigmadphi,2) == size(Ym,2))
                    ddddSdphidphidphidphi = ddddsigmadphidphidphidphi;
                end
            elseif(size(dsigmadphi,2) == size(Ym,2))
                if(size(dsigmadphi,1) == 1)
                    ddddSdphidphidphidphi = repmat(ddddsigmadphidphidphidphi,[size(Ym,1),1,1,1]);
                end
            elseif(and(size(dsigmadphi,1)==1,size(dsigmadphi,2)==1))
                ddddSdphidphidphidphi = repmat(ddddsigmadphidphidphidphi,[size(Ym),1,1,1]);
            end
        end
        
        for k = 1:length(phi) % first order derivatives
            temp = dSdphi(:,:,k);
            dSigmadphi(:,k) = temp(ind);
            if nargout >= 3 % second order derivatives
                for l = 1:length(phi)
                    temp = ddSdphidphi(:,:,k,l);
                    ddSigmadphidphi(:,k,l) = temp(ind);
                    if nargout >= 9 % third order derivatives
                        for m = 1:length(phi)
                            temp = dddSdphidphidphi(:,:,k,l,m);
                            dddSigmadphidphidphi(:,k,l,m) = temp(ind);
                            if nargout >= 14 % fourth order derivatives
                                for n = 1:length(phi)
                                    temp = ddddSdphidphidphidphi(:,:,k,l,m,n);
                                    ddddSigmadphidphidphidphi(:,k,l,m,n) = temp(ind);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % noise model
    switch(Model.exp{s}.noise_model)
        case 'normal'
            if nargout <= 1
                J_D = normal_noise(Y,Ym,Sigma,ind);
            elseif nargout <=2 % first order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma] = normal_noise(Y,Ym,Sigma,ind);
            elseif nargout <=8 % second order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = normal_noise(Y,Ym,Sigma,ind);
            elseif nargout <=13 % third order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                    dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma] = normal_noise(Y,Ym,Sigma,ind);
            else % fourth order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                    dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma,...
                    ddddJ_DdYdYdYdY,ddddJ_DdYdYdYdSigma,ddddJ_DdYdYdSigmadSigma,ddddJ_DdYdSigmadSigmadSigma,ddddJ_DdSigmadSigmadSigmadSigma] = normal_noise(Y,Ym,Sigma,ind);
            end
            
        case 'lognormal'
            
            if nargout <=1
                J_D = lognormal_noise(Y,Ym,Sigma,ind);
            elseif nargout <=2 % first order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma] = lognormal_noise(Y,Ym,Sigma,ind);
            elseif nargout <=8 % second order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = lognormal_noise(Y,Ym,Sigma,ind);
            elseif nargout <=13 % third order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                    dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma] = lognormal_noise(Y,Ym,Sigma,ind);
            else % fourth order derivatives
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                    dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma,...
                    ddddJ_DdYdYdYdY,ddddJ_DdYdYdYdSigma,ddddJ_DdYdYdSigmadSigma,ddddJ_DdYdSigmadSigmadSigma,ddddJ_DdSigmadSigmadSigmadSigma] = lognormal_noise(Y,Ym,Sigma,ind);
            end
            
    end
    
    % parameter model
    switch(Model.exp{s}.parameter_model)
        case 'normal'
            if nargout <=1
                J_b = normal_param(b,delta,type_D);
            elseif nargout <=2 % first order derivatives
                [J_b,dJ_bdb,~]= normal_param(b,delta,type_D);
            elseif nargout <=8 % second order derivatives
                [J_b,dJ_bdb,~,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta] = normal_param(b,delta,type_D);
            elseif nargout <=13 % third order derivatives
                [J_b,dJ_bdb,~,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta,dddJ_bdbdbdb,dddJ_bdbdbddelta,dddJ_bdbddeltaddelta] = normal_param(b,delta,type_D);
            else % fourth order derivatives
                [J_b,dJ_bdb,~,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta,dddJ_bdbdbdb,dddJ_bdbdbddelta,dddJ_bdbddeltaddelta,ddddJ_bdbdbdbdb,ddddJ_bdbdbdbddelta,ddddJ_bdbdbddeltaddelta] = normal_param(b,delta,type_D);
            end
        case 'lognormal'
    end
    
    J = J_D + J_b;
    varargout{1} = J;
    
    if nargout >= 2
        %% dJdb
        dphidb = Model.exp{s}.dphidb(beta,b);
        
        dJ_Ddphi = chainrule_dxdy_dydz(dJ_DdY,dYdphi) + chainrule_dxdy_dydz(dJ_DdSigma,dSigmadphi);
        
        dJ_Ddb = chainrule_dxdy_dydz(dJ_Ddphi,dphidb);
        
        dJdb = dJ_Ddb + dJ_bdb;
        
        varargout{2} = dJdb;
        
        if nargout >= 3
            %% ddJdbdb
            
            ddJ_DdphidY = bsxfun(@times,ddJ_DdYdY,permute(dYdphi,[3,1,2])) + bsxfun(@times,ddJ_DdYdSigma,permute(dSigmadphi,[3,1,2]));
            ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dYdphi,[3,1,2])) + bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigmadphi,[3,1,2]));
            ddJ_Ddphidphi = squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(dYdphi,[3,1,4,2])) + bsxfun(@times,ddJ_DdphidSigma,permute(dSigmadphi,[3,1,4,2])),2));
            
            ddphidbdb = Model.exp{s}.ddphidbdb(beta,b);
            
            ddJ_Ddbdphi = transpose(squeeze(sum(bsxfun(@times,ddJ_Ddphidphi,permute(dphidb,[3,1,2])),2)));
            ddJ_Ddbdb = squeeze(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidb,[3,1,2,4])),2)) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbdb);
            
            ddJdbdb = squeeze(ddJ_Ddbdb) + squeeze(ddJ_bdbdb);
            
            varargout{3} = ddJdbdb;
            
            if nargout >= 4
                
                %% ddJdbdbeta
                dphidbeta = Model.exp{s}.dphidbeta(beta,b);
                ddphidbdbeta = Model.exp{s}.ddphidbdbeta(beta,b);
                
                % if size of b == 1 then we have to permute second term,
                if(numel(b) == 1)
                    ddJ_Ddbdbeta = permute(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidbeta,[3,1,2])),2),[1,3,4,2]) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbdbeta);
                else
                    ddJ_Ddbdbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidbeta,[3,1,2])),2)) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbdbeta);
                end
                
                
                ddJdbdbeta = ddJ_Ddbdbeta;
                
                varargout{4} = ddJdbdbeta;
                
                %% ddJdbdelta
                
                ddJdbddelta = ddJ_bdbddelta;
                
                varargout{5} = ddJdbddelta;
                
                %% ddJdbetadbeta
                
                ddphidbetadbeta = Model.exp{s}.ddphidbetadbeta(beta,b);
                
                ddJ_Ddbetadphi = permute(sum(bsxfun(@times,ddJ_Ddphidphi,permute(dphidbeta,[3,1,2])),2),[3,2,1]);
                if(numel(b)==1)
                    
                else
                    ddJ_Ddbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbetadbeta);
                end
                
                ddJdbetadbeta = ddJ_Ddbetadbeta;
                
                varargout{6} = ddJdbetadbeta;
                
                %% ddJddeltaddelta
                
                ddJddeltaddelta = ddJ_bddeltaddelta;
                
                varargout{7} = ddJddeltaddelta;
                
                %% ddJdbetaddelta
                
                ddJdbetaddelta = zeros(length(beta),size(dDddelta,3));
                
                varargout{8} = ddJdbetaddelta;
                
                if nargout >= 9
                    
                    %% dddJdbdbdb
                    
                    temp = squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(ddYdphidphi,[4,1,5,2,3])) ...
                        + bsxfun(@times,ddJ_DdphidSigma,permute(ddSigmadphidphi,[4,1,5,2,3])),2));
                    
                    dddJ_Ddphidphidphi = chainrule_dxdy_dydz(dJ_DdSigma,dddSigmadphidphidphi) ...
                        + temp + permute(temp,[2,1,3]) + permute(temp,[2,3,1]);
                    
                    
                    dddJ_DdphidYdY = bsxfun(@times,dddJ_DdYdYdY,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_DdYdYdSigma,permute(dSigmadphi,[3,1,2]));
                    dddJ_DdphidYdSigma = bsxfun(@times,dddJ_DdYdYdSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_DdYdSigmadSigma,permute(dSigmadphi,[3,1,2]));
                    dddJ_DdphidSigmadSigma = bsxfun(@times,dddJ_DdYdSigmadSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_DdSigmadSigmadSigma,permute(dSigmadphi,[3,1,2]));
                    dddJ_DdphidphidY = bsxfun(@times,dddJ_DdphidYdY,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_DdphidYdSigma,permute(dSigmadphi,[3,1,4,2]));
                    dddJ_DdphidphidSigma = bsxfun(@times,dddJ_DdphidYdSigma,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_DdphidSigmadSigma,permute(dSigmadphi,[3,1,4,2]));
                    
                    dddJ_Ddphidphidphi = dddJ_Ddphidphidphi + squeeze(sum(bsxfun(@times,dddJ_DdphidphidY,permute(dYdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,dddJ_DdphidphidSigma,permute(dSigmadphi,[3,1,4,5,2])),2));
                    
                    dddJ_Ddbdphidphi = permute(sum(bsxfun(@times,dddJ_Ddphidphidphi,permute(dphidb,[1,3,4,2])),1),[4,2,3,1]);
                    dddJ_Ddbdbdphi = permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(dphidb,[3,1,4,2])),2),[1,4,3,2]) ...
                        + permute(sum(bsxfun(@times,ddJ_Ddphidphi,permute(ddphidbdb,[1,4,2,3])),1),[3,4,2,1]);
                    dddJ_Ddbdbdb = permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(dphidb,[3,4,1,2])),3),[1,2,4,3]) ...
                        + squeeze(chainrule_dxdy_dydz(ddJ_Ddbdphi,ddphidbdb));
                    
                    dddJdbdbdb = dddJ_Ddbdbdb + squeeze(dddJ_bdbdbdb);
                    
                    varargout{9} = dddJdbdbdb;
                    
                    %% dddJdbdbbeta
                    
                    dddJ_Ddbdbdbeta = permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                    
                    dddJdbdbdbeta = dddJ_Ddbdbdbeta;
                    
                    varargout{10} = dddJdbdbdbeta;
                    
                    
                    %% dddJdbdbdelta
                    dddJdbdbddelta = dddJ_bdbdbddelta;
                    
                    varargout{11} = squeeze(dddJdbdbddelta);
                    
                    %% dddJdbdbetadbeta
                    
                    dddJ_Ddbdbetadphi = permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(dphidbeta,[3,1,4,2])),2),[1,4,3,2]);
                    dddJ_Ddbdbetadbeta = permute(sum(bsxfun(@times,dddJ_Ddbdbetadphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                    
                    dddJdbdbetadbeta = dddJ_Ddbdbetadbeta;
                    
                    varargout{12} = squeeze(dddJdbdbetadbeta);
                    
                    %% dddJdbdbdelta
                    dddJdbddeltaddelta = dddJ_bdbddeltaddelta;
                    
                    varargout{13} = squeeze(dddJdbddeltaddelta);
                    
                    
                    if nargout >= 14
                        
                        %% ddddJdbdbdbdb
                        
                        temp = squeeze(sum(bsxfun(@times,bsxfun(@times,ddJ_DdYdY,permute(ddYdphidphi,[4,1,2,3])),permute(ddYdphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_DdYdSigma,permute(ddYdphidphi,[4,1,2,3])),permute(ddSigmadphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_DdSigmadSigma,permute(ddSigmadphidphi,[4,1,2,3])),permute(ddSigmadphidphi,[4,1,5,6,2,3])),2));
                        
                        ddddJ_Ddphidphidphidphi = ...
                            permute(temp,[1,3,2,4]) + permute(temp,[1,3,4,2]) + permute(temp,[3,4,1,2]);
                        
                        ddddJ_DdphidYdYdY = bsxfun(@times,ddddJ_DdYdYdYdY,permute(dYdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_DdYdYdYdSigma,permute(dSigmadphi,[3,1,2]));
                        ddddJ_DdphidYdYdSigma = bsxfun(@times,ddddJ_DdYdYdYdSigma,permute(dYdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_DdYdYdSigmadSigma,permute(dSigmadphi,[3,1,2]));
                        ddddJ_DdphidYdSigmadSigma = bsxfun(@times,ddddJ_DdYdYdSigmadSigma,permute(dYdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_DdYdSigmadSigmadSigma,permute(dSigmadphi,[3,1,2]));
                        ddddJ_DdphidSigmadSigmadSigma = bsxfun(@times,ddddJ_DdYdSigmadSigmadSigma,permute(dYdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_DdSigmadSigmadSigmadSigma,permute(dSigmadphi,[3,1,2]));
                        ddddJ_DdphidphidYdY = bsxfun(@times,ddddJ_DdphidYdYdY,permute(dYdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_DdphidYdYdSigma,permute(dSigmadphi,[3,1,4,2]));
                        ddddJ_DdphidphidYdSigma = bsxfun(@times,ddddJ_DdphidYdYdSigma,permute(dYdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_DdphidYdSigmadSigma,permute(dSigmadphi,[3,1,4,2]));
                        ddddJ_DdphidphidSigmadSigma = bsxfun(@times,ddddJ_DdphidYdSigmadSigma,permute(dYdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_DdphidSigmadSigmadSigma,permute(dSigmadphi,[3,1,4,2]));
                        ddddJ_DdphidphidphidY = bsxfun(@times,ddddJ_DdphidphidYdY,permute(dYdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_DdphidphidYdSigma,permute(dSigmadphi,[3,1,4,5,2]));
                        ddddJ_DdphidphidphidSigma = bsxfun(@times,ddddJ_DdphidphidYdSigma,permute(dYdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_DdphidphidSigmadSigma,permute(dSigmadphi,[3,1,4,5,2]));
                        ddddJ_Ddphidphidphidphi = ddddJ_Ddphidphidphidphi + squeeze(sum(bsxfun(@times,ddddJ_DdphidphidphidY,permute(dYdphi,[3,1,4,5,6,2])) ...
                            + bsxfun(@times,ddddJ_DdphidphidphidSigma,permute(dSigmadphi,[3,1,4,5,6,2])),2));
                        
                        ddddJ_Ddbdphidphidphi = permute(sum(bsxfun(@times,ddddJ_Ddphidphidphidphi,permute(dphidb,[1,3,4,5,2])),1),[5,2,3,4,1]);
                        ddddJ_Ddbdbdphidphi = permute(sum(bsxfun(@times,ddddJ_Ddbdphidphidphi,permute(dphidb,[3,1,4,5,2])),2),[1,5,3,4,2]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddphidphidphi,permute(ddphidbdb,[1,4,5,2,3])),1),[4,5,2,3,1]);
                        ddddJ_Ddbdbdbdphi = permute(sum(bsxfun(@times,ddddJ_Ddbdbdphidphi,permute(dphidb,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(ddphidbdb,[4,1,5,2,3])),2),[1,4,5,3,2]);
                        ddddJ_Ddbdbdbdb = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbdphi,permute(dphidb,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbdb,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        ddddJdbdbdbdb = ddddJ_Ddbdbdbdb + squeeze(ddddJ_bdbdbdbdb);
                        
                        varargout{14} = ddddJdbdbdbdb;
                        
                        %% ddddJdbdbdbdbeta
                        
                        ddddJ_Ddbdbdbdbeta = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbdphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        ddddJdbdbdbdbeta = ddddJ_Ddbdbdbdbeta;
                        
                        varargout{15} = ddddJdbdbdbdbeta;
                        
                        %% ddddJdbdbdbetadbeta
                        
                        ddddJ_Ddbdbdbetadphi = permute(sum(bsxfun(@times,ddddJ_Ddbdbdphidphi,permute(dphidbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(ddphidbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
                        ddddJ_Ddbdbdbetadbeta = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbetadphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        ddddJdbdbdbetadbeta = ddddJ_Ddbdbdbetadbeta;
                        
                        varargout{16} = ddddJdbdbdbetadbeta;
                        
                        %% ddddJdbdbdbddelta
                        
                        ddddJdbdbdbddelta = ddddJ_bdbdbdbddelta;
                        
                        varargout{17} = ddddJdbdbdbddelta;
                        
                        %% ddddJdbdbddeltaddelta
                        
                        ddddJdbdbddeltaddelta = ddddJ_bdbdbddeltaddelta;
                        
                        varargout{18} = ddddJdbdbddeltaddelta;
                        
                        %% dddddJdbdbdbetaddelta
                        varargout{19} = zeros(length(b),length(b),length(beta),length(b));
                        
                        
                    end
                end
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR SIGMA POINTS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,sy] = simulateForSP(model,tout,phi,kappa)
    
    % Simulate model
    [~,~,~,y,~,sy] = model(tout,phi,kappa);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_noise(Y,Ym,Sigma,ind)
    if nargout >=1
        % J_D
        varargout{1} = sum(0.5*((Y(ind) - Ym(ind))./Sigma(ind)).^2 + 0.5*log(sqrt(2*pi)*Sigma(ind).^2));
        if nargout >= 3
            % dJ_DdY
            varargout{2} = transpose((Y(ind) - Ym(ind))./(Sigma(ind).^2));
            % dJ_DdSigma
            varargout{3} = transpose(- (((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^3)) + 1./Sigma(ind));
            if nargout >= 4
                %ddJ_DdYdY
                varargout{4} = transpose(1./(Sigma(ind).^2));
                %ddJ_DdYdSigma
                varargout{5} = transpose(-2*(Y(ind) - Ym(ind))./(Sigma(ind).^3));
                %ddJ_DdSigmadSigma
                varargout{6} = transpose(3*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^4)) - 1./(Sigma(ind).^2));
                if nargout >= 7
                    %dddJ_DdYdYdY
                    varargout{7} = transpose(zeros(size(Y(ind))));
                    %dddJ_DdYdYdSigma
                    varargout{8} = transpose(- 2./(Sigma(ind).^3));
                    %dddJ_DdYdSigmadSigma
                    varargout{9} = transpose(6*(Y(ind) - Ym(ind))./(Sigma(ind).^4));
                    %dddJ_DdSigmadSigmadSigma
                    varargout{10} = transpose(- 12*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^5)) + 2./(Sigma(ind).^3));
                    if nargout >= 11
                        %ddddJ_DdYdYdYdY
                        varargout{11} = transpose(zeros(size(Y(ind))));
                        %ddddJ_DdYdYdYdSigma
                        varargout{12} = transpose(zeros(size(Y(ind))));
                        %ddddJ_DdYdYdSigmadSigma
                        varargout{13} = transpose(6./(Sigma(ind).^4));
                        %ddddJ_DdYdSigmadSigmadSigma
                        varargout{14} = transpose(- 24*((Y(ind) - Ym(ind))./(Sigma(ind).^5)));
                        %ddddJ_DdSigmadSigmadSigmadSigma
                        varargout{15} = transpose(60*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^5)) - 6./(Sigma(ind).^4));
                    end
                end
            end
        end
    end
end

function varargout = lognormal_noise(Y,Ym,Sigma)
    if nargout >=1
        % J_D
        varargout{1} = transpose(sum(0.5*((log(Y) - Ym)./Sigma).^2 + 0.5*log(sqrt(2*pi)*(Sigma.^2).*(Y.^2))));
        if nargout >= 3
            % dJ_DdY
            varargout{2} = transpose((log(Y) - Ym)./(Sigma.^2).*(1./Y) + 1./Y);
            % dJ_DdSigma
            varargout{3} = transpose(- (((log(Y) - Ym).^2)./(Sigma.^3)) + 1./Sigma);
            if nargout >= 4
                %ddJ_DdYdY
                varargout{4} = transpose((1-(log(Y) - Ym))./(Sigma.^2.*Y.^2) - 1./(Y.^2));
                %ddJ_DdYdSigma
                varargout{5} = transpose(-2*(Y - Ym)./(Sigma.^3.*Y));
                %ddJ_DdSigmadSigma
                varargout{6} = transpose(3*(((Y - Ym).^2)./(Sigma.^4)) - 1./(Sigma^2));
                if nargout >= 7
                    %dddJ_DdYdYdY
                    varargout{7} = transpose(- 3./(Sigma.^2.*Y.^3)-2./(Y.^3));
                    %dddJ_DdYdYdSigma
                    varargout{8} = transpose((1-(log(Y) - Ym))./(Sigma.^3.*Y.^2));
                    %dddJ_DdYdSigmadSigma
                    varargout{9} = transpose(+ 6*(Y - Ym)./(Sigma.^4.*Y));
                    %dddJ_DdSigmadSigmadSigma
                    varargout{10} = transpose(- 12*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^5)) + 2./(Sigma(ind)^3));
                    if nargout >= 11
                        %ddddJ_DdYdYdYdY
                        varargout{11} = transpose(+9./(Sigma.^2.*Y.^4)+6./(Y.^4));
                        %ddddJ_DdYdYdYdSigma
                        varargout{12} = transpose(6./(Sigma.^3.*Y.^3)-2./(Y.^3));
                        %ddddJ_DdYdYdSigmadSigma
                        varargout{13} = transpose(-3*(1-(log(Y) - Ym))./(Sigma.^4.*Y.^2));
                        %ddddJ_DdYdSigmadSigmadSigma
                        varargout{14} = transpose(- 24*(Y - Ym)./(Sigma.^5.*Y));
                        %ddddJ_DdSigmadSigmadSigmadSigma
                        varargout{15} = transpose(- 60*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^6)) + 6./(Sigma(ind)^4));
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR PARAMETER DENSITIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_param(b,delta,type_D)

    [D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);

    if nargout >=1
        % J_b
        varargout{1} = 0.5*b'*invD*b ...
            + 0.5*log(det(D));
        if nargout >= 2
            % dJ_bdb
            varargout{2} = transpose(invD*b);
            % dJ_bddelta
            varargout{3} = transpose(0.5*squeeze(sum(sum(bsxfun(@times,dinvDddelta,bsxfun(@times,permute(b,[2,1]),permute(b,[1,2]))),1),2)) ... % 1/2*b*dinvD*b
                +0.5*squeeze(sum(sum(sum(bsxfun(@times,invD.*eye(length(b)),permute(dDddelta,[4,1,2,3])),2),1),3))); % 1/2*Tr(invD*dD)
            if nargout >= 4
                % ddJ_bdbdb
                varargout{4} = invD;
                % ddJ_bdbddelta
                varargout{5} = squeeze(sum(bsxfun(@times,dinvDddelta,b),2)); % dinvD*b it does not reall matter onto which of the first two dimensions we multiply b here as D and all its derivatives are symmetric
                % ddJ_bddeltaddelta
                varargout{6} = 0.5*squeeze(sum(sum(bsxfun(@times,ddinvDddeltaddelta,bsxfun(@times,permute(b,[2,1]),permute(b,[1,2]))),1),2)) ... % 1/2*b*ddinvD*b
                    +0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2)),eye(length(b))),1),2)); % 1/2*Tr(dinvD*dD + invD*ddD)
                if nargout >= 7
                    % dddJ_bdbdbdb
                    varargout{7} = zeros(length(b),length(b),length(b));
                    % dddJ_bdbdbddelta
                    varargout{8} = dinvDddelta;
                    % dddJ_bdbddeltaddelta
                    varargout{9} = squeeze(sum(bsxfun(@times,ddinvDddeltaddelta,b),2)); % dinvD*b it does not reall matter onto which of the first two dimensions we multiply b here as D and all its derivatives are symmetric
                    if nargout >= 10
                        % ddddJ_bdbdbdbdb
                        varargout{10} = zeros(length(b),length(b),length(b),length(b));
                        % ddddJ_bdbdbdbddelta
                        varargout{11} = zeros(length(b),length(b),length(b),length(b));
                        % ddddJ_bdbdbddeltaddelta
                        varargout{12} = ddinvDddeltaddelta;
                        
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR PARAMETER PENALTY   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = penal_param(b_s,delta,type_D)


    [D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);
    
    
    d_s = size(b_s,2);
    n_b = size(b_s,1);
    
    mu_s = 1/d_s*(sum(b_s,2));
    
    S_s = (b_s*b_s') + 1e-10*eye(n_b);
    
    
    
    if nargout >=1
        % J_s = -log(p({mu_S,S_s}|D))
        varargout{1} = + 0.5*(d_s-n_b-1)*log(det(S_s)) ... % log(det(S))
            - 0.5*(d_s+1)*log(det(D)) ... % log(det(D))
            - 0.5*d_s*(mu_s'*D*mu_s) ... % mu*D*mu
            - 0.5*trace(invD*S_s); % tr(invD*S)
        
        if nargout >= 2
            dmu_sdb_s = 1/d_s*repmat(eye(n_b),[1,1,d_s]);
            dS_sdb_s = bsxfun(@times,permute(b_s,[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) + bsxfun(@times,permute(b_s,[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
            invS_s = pinv(S_s);
            
            % dJ_sdb_s
            varargout{2} = 0.5*(d_s-n_b-1)*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invS_s,permute(dS_sdb_s,[5,1,2,3,4])),2)),1),3)) ... % tr(invS*dS)
                - d_s*squeeze(sum(bsxfun(@times,mu_s'*D,permute(dmu_sdb_s,[4,1,2,3])),2)) ... % mu*D*dmu
                - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dS_sdb_s,[5,1,2,3,4])),2)),1),3)); % tr(invD*dS)
            
            % dJ_sddelta
            varargout{3} = transpose(-0.5*(d_s+1)*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dDddelta,[4,1,2,3])),2)),1),3)) ... % tr(invD*dD)
                - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,dDddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)) ... % mu*dD*mu
                - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(S_s,[3,1,2]),permute(dinvDddelta,[1,2,4,3])),2)),1),3))); % tr(dinvD*S)
            if nargout >= 4
                ddS_sdb_sdb_s = bsxfun(@times,permute(ones(size(b_s)),[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) ...
                    + bsxfun(@times,permute(ones(size(b_s)),[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
                dinvS_sdb_s = -squeeze(sum(sum(bsxfun(@times,bsxfun(@times,invS_s,permute(dS_sdb_s,[5,1,2,6,3,4])),permute(invS_s,[3,4,1,2])),2),3));
                % ddJ_sdb_sdb_s
                varargout{4} = 0.5*(d_s-n_b-1)*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvS_sdb_s,[1,2,5,3,4]),permute(dS_sdb_s,[5,1,2,3,4])),2)+sum(bsxfun(@times,invS_s,permute(ddS_sdb_sdb_s,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ... % 1/2*Tr(dinvS*dS + invS*ddS)
                    - d_s*squeeze(sum(sum(bsxfun(@times,bsxfun(@times,permute(dmu_sdb_s,[1,4,2,3]),D),permute(dmu_sdb_s,[4,1,2,3])),2),1)) ...  % dmu*D*dmu
                    - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(ddS_sdb_sdb_s,[5,1,2,3,4])),2)),1),3)); % tr(invD*ddS)
                
                
                % ddJ_sdb_sddelta
                varargout{5} = - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,permute(dDddelta,[1,2,4,5,3]),bsxfun(@times,permute(mu_s,[2,1]),permute(dmu_sdb_s,[1,4,2,3]))),1),2)) ... % mu*dD*dmu
                    - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(dS_sdb_s,[1,2,5,3,4]),permute(dinvDddelta,[4,1,2,5,6,3])),2)),3),1)); % tr(dinvD*dS)
                
                % ddJ_sddeltaddelta
                varargout{6} = -0.5*(d_s+1)*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ...% 1/2*Tr(ddinvD*dD + invD*ddD)
                    - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,ddDddeltaddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)) ... % mu*ddD*mu
                    - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(S_s,[3,1,2]),permute(ddinvDddeltaddelta,[1,2,5,3,4])),2)),1),3)); % tr(ddinvD*S)
                
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR COMPUTATION OF DERIVATIVES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdz = chainrule_dxdy_dydz(dxdy,dydz)
    d1 = ndims(dxdy);
    d2 = ndims(dydz);
    
    if(~size(dxdy,d1)==size(dydz,1))
        error('dimensions must agree')
    end
    %      dx        dy
    % [ 1 : d1 - 1 , d1     , d1 + 1:d2-1]
    %               <.,.>
    %   ********     dy          dz
    % [ 1 : d1 - 1 , d1     , d1 + 1:d2-1]
    if(numel(dxdy)>1)
        if(d1>1)
            if(and(d1==2,d2==2))
                dxdz = dxdy*dydz;
            else
                dxdz = permute(sum(bsxfun(@times,dxdy,permute(dydz,[d2+(1:(d1-1)),1,2:d2])),d1),[1:(d1-1),(d1+1):(d1+d2-1),d1]);
            end
        else
            dxdz =     permute(sum(bsxfun(@times,dxdy,permute(dydz,[              1,2:d2])),d1),[1:(d1-1),(d1+1):(d1+d2-1),d1]);
        end
    else
        dxdz = permute(dxdy*dydz,[d2+(1:(d1-1)),2:d2,1]);
    end    
    dxdz = squeeze(dxdz);
end

function ddxdzdz = chainrule_ddxdydy_dydz(ddxdydy,dydz)
    d1 = ndims(ddxdydy);
    d2 = ndims(dydz);
    %      dx        dy      dy
    % [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d2-1]
    %               <.,.>
    %   ********     dy     ***           dz       **************
    % [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d2-1]
    %   ******************   dy      **********         dz
    %                       <.,.>
    % [ 1 : d1 - 1        , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d2-1]
    if(d1>2)
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[d2+(1:(d1-2)),1,d2+d1-1,2:d2             ]     )),d1-1);
    else
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[              1,d2+d1-1,2:d2             ]     )),d1-1);
    end
    ddxdzdz =     squeeze(sum(bsxfun(@times,ddxdzdy,permute(dydz,[d2+(1:(d1-1))  ,1      ,d2+(d1:(d1+d2-1)),2:d2])),d1));
    
end

function ddxdzdv = chainrule_ddxdydy_dydz_dydv(ddxdydy,dydz,dydv)
    d1 = ndims(ddxdydy);
    d2 = ndims(dydz);
    d3 = ndims(dydv);
    %      dx        dy      dy
    % [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d3-1]
    %               <.,.>
    %   ********     dy     ***           dz       **************
    % [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d3-1]
    %   ******************   dy      **********         dv
    %                       <.,.>
    % [ 1 : d1 - 1        , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d3-1]
    if(d1>2)
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[d2+(1:(d1-2)),1,d2+d1-1,2:d2             ]     )),d1-1);
    else
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[              1,d2+d1-1,2:d2             ]     )),d1-1);
    end
    ddxdzdv =     squeeze(sum(bsxfun(@times,ddxdzdy,permute(dydv,[d3+(1:(d1-1))  ,1      ,d3+(d1:(d1+d2-1)),2:d3])),d1));
    
end

% function dddxdzdzdz = chainrule_dddxdydydy_dydz(dddxdydydy,dydz)
%     d1 = ndims(dddxdydydy);
%     d2 = ndims(dydz);
%     %      dx        dy      dy       dy
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d2-1 ]
%     %               <.,.>
%     %   ********     dy    *************           dz
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d2-1 ]
%     %                        <.,.>
%     %   ******************   dy    *********************         dz
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d2-1 ]
%         %                            <.,.>
%     %   ***************************   dy    ********************************         dz
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d2-1 ]
%     if(d1>1)
%         dddxdzdydy =         sum(bsxfun(@times,dddxdydydy,permute(dydz,[d2+(1:(d1-3)),1,d2+(d1-2:d1-1) ,2:d2 ]                   )),d1-2);
%     else
%         dddxdzdydy =         sum(bsxfun(@times,dddxdydydy,permute(dydz,[              1,d2+(d1-2:d1-1) ,2:d2 ]                   )),d1-2);
%     end
%     dddxdzdzdy =             sum(bsxfun(@times,dddxdzdydy,permute(dydz,[d2+(1:(d1-2))  ,1,d2+(d1-1:d2+(d1-1)),2:d2         ]     )),d1-1);
%     dddxdzdzdz =     squeeze(sum(bsxfun(@times,dddxdzdzdy,permute(dydz,[d2+(1:(d1-1))    ,1            ,d2+(d1:2*d2+(d1-1)),2:d2])),d1));
%
% end

% function dddxdzdzdv = chainrule_dddxdydydv_dydz_dydv(dddxdydydy,dydz,dydv)
%     d1 = ndims(dddxdydydy);
%     d2 = ndims(dydz);
%     d3 = ndims(dydv);
%     %      dx        dy      dy       dy
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d3-1 ]
%     %               <.,.>
%     %   ********     dy    *************           dz
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d3-1 ]
%     %                        <.,.>
%     %   ******************   dy    *********************         dz
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d3-1 ]
%         %                            <.,.>
%     %   ***************************   dy    ********************************         dv
%     % [ 1 : d1-3 , d1 - 2 , d1 - 1 ,  d1 ,  d1 + 1 : d2-1 , d1 + d2 + 1:d2-1, d1 + 2*d2 + 1:d3-1 ]
%     if(d1>1)
%         dddxdzdydy =     sum(bsxfun(@times,dddxdydydy,permute(dydz,[d2+(1:(d1-3)),1,d2+(d1-2:d1-1) ,2:d2 ]                   )),d1-2);
%     else
%         dddxdzdydy =     sum(bsxfun(@times,dddxdydydy,permute(dydz,[              1,d2+(d1-2:d1-1) ,2:d2 ]                   )),d1-2);
%     end
%     dddxdzdzdy =         sum(bsxfun(@times,dddxdzdydy,permute(dydz,[d2+(1:(d1-2))  ,1,d2+(d1-1:d2+(d1-1)),2:d2         ]     )),d1-1);
%     dddxdzdzdv = squeeze(sum(bsxfun(@times,dddxdzdzdy,permute(dydv,[d3+(1:(d1-1))    ,1            ,d3+(d1:2*d2+(d1-1)),2:d3])),d1));
%
% end

% function ddddxdzdzdzdz = chainrule_ddddxdydydydy_dydz(ddddxdydydydy,dydz)
%     d1 = ndims(ddddxdydydydy);
%     d2 = ndims(dydz);
%
%     if(d1>4)
%         ddddxdzdydydy =     sum(bsxfun(@times,ddddxdydydydy,permute(dydz,[d2+(1:(d1-4)),1,d2+((d1-3):(d1-1)),2:d2]                      )),d1-3);
%     else
%         ddddxdzdydydy =     sum(bsxfun(@times,ddddxdydydydy,permute(dydz,[              1,d2+((d1-3):(d1-1)),2:d2]                      )),d1-3);
%     end
%     ddddxdzdzdydy =         sum(bsxfun(@times,ddddxdzdydydy,permute(dydz,[d2+(1:(d1-3))  ,1,d2+((d1-2):d2+(d1-1)),2:d2]                 )),d1-2);
%     ddddxdzdzdzdy =         sum(bsxfun(@times,ddddxdzdzdydy,permute(dydz,[d2+(1:(d1-2))    ,1,d2+((d1-1):2*d2+(d1-1)) ,2:d2    ]        )),d1-1);
%     ddddxdzdzdzdz = squeeze(sum(bsxfun(@times,ddddxdzdzdzdy,permute(dydz,[d2+(1:(d1-1))      ,1             ,d2+(d1:3*d2+(d1-1)),2:d2-1])),d1));
%
% end




