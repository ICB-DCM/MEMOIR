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
        logL_T = 0;
        logL_b = 0;
        logL_I = 0;
        
        % Evaluation of time index set
        [~,ind_time] = ismember(Data{s}.SCTL.time,t_s);
        
        % Initialization
        Sim_SCTL.Y = nan(size(Data{s}.SCTL.Y));
        Sim_SCTL.T = nan(size(Data{s}.SCTL.T));
        Sim_SCTL.R = nan(size(Data{s}.SCTL.T));
        
        % Loop: Indiviudal cells
        for i = 1:size(Data{s}.SCTL.Y,3)
            % Load single-cell data
            Ym_si = Data{s}.SCTL.Y(ind_time,:,i);
            Tm_si = Data{s}.SCTL.T(:,:,i);
            
            ind_y = find(~isnan(Ym_si));
            ind_t = find(~isnan(Tm_si));
            
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
                        = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
                case 1
                    F_diff = 0;
                    b_diff = 0;
                    [bhat_si,G] ...
                        = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
                case 2
                    b_diff = 1;
                    if(Model.integration)
                        F_diff = 3;
                        
                        [bhat_si,dbhat_sidbeta,dbhat_siddelta,...
                            G,dGdb,pdGpdbeta,pdGpddelta] ...
                            = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,2)
                    else
                        F_diff = 2;
                        [bhat_si,dbhat_sidbeta,dbhat_siddelta,...
                            G] ...
                            = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
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
                            = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
                    else
                        F_diff = 3;
                        b_diff = 2;
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,2)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,3)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,2,4)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,2,5)
                        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,3,6)
                        [bhat_si,dbhat_sidbeta,dbhat_siddelta,ddbhat_sidbetadbeta,ddbhat_sidbetaddelta,ddbhat_siddeltaddelta,...
                            G,dGdb,pdGpdbeta] ...
                            = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
                    end
            end
            
            % Store bhat
            P{s}.SCTL.bhat(:,i) = bhat_si;
            
            % Construct single-cell parameter
            phi_si = Model.exp{s}.phi(beta,bhat_si);
            
            % Simulate model
            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(phi_si,@(phi) Model.exp{s}.model(t_s,phi,Data{s}.condition),1e-4,4,6)
            if nargout <2
                option_simu.sensi = 0;
                sol = Model.exp{s}.model(t_s,phi_si,Data{s}.condition,option_simu);
            elseif nargout <3
                option_simu.sensi = 1;
                sol = Model.exp{s}.model(t_s,phi_si,Data{s}.condition,option_simu);
                dYdphi = sol.sy;
                dTdphi = sol.sroot(ind_t,:,:);
                dRdphi = sol.srootval(ind_t,:,:);
            else
                if(Model.integration)
                    option_simu.sensi = 1;
                    sol = Model.exp{s}.model(t_s,phi_si,Data{s}.condition,option_simu);
                    dYdphi = sol.sy;
                    dTdphi = sol.sroot(ind_t,:,:);
                    dRdphi = sol.srootval(ind_t,:,:);
                else
                    option_simu.sensi = 1;
                    sol = Model.exp{s}.model(t_s,phi_si,Data{s}.condition,option_simu);
                    dYdphi = sol.sy;
                    dTdphi = sol.sroot(ind_t,:,:);
                    dRdphi = sol.srootval(ind_t,:,:);
                end
            end
            Y = sol.y;
            T = sol.root;
            R = sol.rootval;
            Y_si = Y(ind_time,:);
            T_si = T(ind_t,:);
            R_si = R(ind_t,:);
            
            
            % Apply indexing to derivatives of Y
            if nargout >= 2
                dY_sidphi = zeros(length(ind_y),length(phi_si));
                dT_sidphi = zeros(length(ind_t),length(phi_si));
                if nargout >= 3
                    ddY_sidphidphi = zeros(length(ind_y),length(phi_si),length(phi_si));
                    ddT_sidphidphi = zeros(length(ind_t),length(phi_si),length(phi_si));
                end
                for k = 1:length(phi_si)
                    tempy = dYdphi(:,:,k);
                    tempt = dTdphi(:,:,k);
                    dY_sidphi(:,k) = tempy(ind_y);
                    dT_sidphi(:,k) = tempt(ind_t);
                    if Model.integration
                        for j = 1:length(phi_si)
                            tempy = ddYdphidphi(:,:,j,k);
                            tempt = ddTdphidphi(:,:,j,k);
                            ddY_sidphidphi(:,j,k) = tempy(ind_y);
                            ddY_sidphidphi(:,j,k) = tempt(ind_t);
                        end
                    end
                end
            end
            
            % Construct sigma
            sigma_noise = Model.exp{s}.sigma_noise(phi_si);
            
            if(size(sigma_noise,1) == size(Ym_si,1))
                if(size(sigma_noise,2) == 1)
                    Sigma_noise_si = repmat(sigma_noise,[1,size(Ym_si,2)]);
                elseif(size(sigma,2) == size(Ym_si,2))
                    Sigma_noise_si = sigma_noise;
                else
                    error('Incompatible size of sigma_noise parametrisation!')
                end
            elseif(size(sigma_noise,2) == size(Ym_si,2))
                if(size(sigma_noise,1) == 1)
                    Sigma_noise_si = repmat(sigma_noise,[size(Ym_si,1),1]);
                else
                    error('Incompatible size of sigma_noise parametrisation!')
                end
            elseif(and(size(sigma_noise,1)==1,size(sigma_noise,2)==1))
                Sigma_noise_si = repmat(sigma_noise,size(Ym_si));
            else
                error('Incompatible size of sigma_noise parametrisation!')
            end
            
            sigma_time = Model.exp{s}.sigma_time(phi_si);
            
            if(size(sigma_time,1) == size(Tm_si,1))
                if(size(sigma_time,2) == 1)
                    Sigma_time_si = repmat(sigma_time,[1,size(Tm_si,2)]);
                elseif(size(sigma,2) == size(Tm_si,2))
                    Sigma_time_si = sigma_time;
                else
                    error('Incompatible size of sigma_time parametrisation!')
                end
            elseif(size(sigma_time,2) == size(Tm_si,2))
                if(size(sigma_time,1) == 1)
                    Sigma_time_si = repmat(sigma_time,[size(Tm_si,1),1]);
                else
                    error('Incompatible size of sigma_time parametrisation!')
                end
            elseif(and(size(sigma_time,1)==1,size(sigma_time,2)==1))
                Sigma_time_si = repmat(sigma_time,size(Tm_si));
            else
                error('Incompatible size of sigma_time parametrisation!')
            end
            
            
            % Evaluation of likelihood and likelihood gradient
            
            % this is part accounts for the noise model
            % J_D = log(p(Y(b,beta)|D))
            switch(Model.exp{s}.noise_model)
                case 'normal'
                    switch(nargout)
                        case 0
                            J_D = normal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                        case 1
                            J_D = normal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                        case 2
                            [J_D,dJ_DdY,dJ_DdSigma] = normal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                        case 3
                            [J_D,dJ_DdY,dJ_DdSigma,ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = normal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                    end
                    
                case 'lognormal'
                    switch(nargout)
                        case 0
                            J_D = lognormal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                        case 1
                            J_D = lognormal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                        case 2
                            [J_D,dJ_DdY,dJ_DdSigma] = lognormal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                        case 3
                            [J_D,dJ_DdY,dJ_DdSigma,ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = lognormal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                    end
            end
            
            % this is part accounts for the time model
            % J_D = log(p(Y(b,beta)|D))
            switch(Model.exp{s}.noise_model)
                case 'normal'
                    switch(nargout)
                        case 0
                            J_T = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                        case 1
                            J_T = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                        case 2
                            [J_T,dJ_TdT,dJ_TdSigma] = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
                        case 3
                            [J_T,dJ_TdT,dJ_TdSigma,ddJ_TdTdT,ddJ_TdTdSigma,ddJ_TdSigmadSigma] = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
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
                        case 3
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(bhat_si,@(b) normal_param(b,delta,type_D),1e-4,1,2)
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,type_D),1e-4,1,3)
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(bhat_si,@(b) normal_param(b,delta,type_D),1e-4,2,4)
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,type_D),1e-4,2,5)
                            % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta) normal_param(bhat_si,delta,type_D),1e-4,3,6)
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
                dsigma_noisedphi = Model.exp{s}.dsigma_noisedphi(phi_si);
                dSigma_noisedphi = zeros(length(ind_y),length(phi_si));
                
                if(size(dsigma_noisedphi,1) == size(Ym_si,1))
                    if(size(dsigma_noisedphi,2) == 1)
                        dSNdphi = repmat(dsigma_noisedphi,[1,size(Ym_si,2),1]);
                    elseif(size(dsigma_noisedphi,2) == size(Ym_si,2))
                        dSNdphi = dsigma_noisedphi;
                    end
                elseif(size(dsigma_noisedphi,2) == size(Ym_si,2))
                    if(size(dsigma_noisedphi,1) == 1)
                        dSNdphi = repmat(dsigma_noisedphi,[size(Ym_si,1),1,1]);
                    end
                elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
                    dSNdphi = repmat(dsigma_noisedphi,[size(Ym_si),1]);
                end
                for k = 1:length(phi_si)
                    temp = dSNdphi(:,:,k);
                    dSigma_noisedphi(:,k) = temp(ind_y);
                end
            end
            
            if nargout >= 2
                dsigma_timedphi = Model.exp{s}.dsigma_timedphi(phi_si);
                dSigma_timedphi = zeros(length(ind_t),length(phi_si));
                
                if(size(dsigma_timedphi,1) == size(Tm_si,1))
                    if(size(dsigma_timedphi,2) == 1)
                        dSTdphi = repmat(dsigma_timedphi,[1,size(Tm_si,2),1]);
                    elseif(size(dsigma_timedphi,2) == size(Tm_si,2))
                        dSTdphi = dsigma_timedphi;
                    end
                elseif(size(dsigma_timedphi,2) == size(Tm_si,2))
                    if(size(dsigma_timedphi,1) == 1)
                        dSTdphi = repmat(dsigma_timedphi,[size(Tm_si,1),1,1]);
                    end
                elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
                    dSTdphi = repmat(dsigma_timedphi,[size(Tm_si),1]);
                end
                for k = 1:length(phi_si)
                    temp = dSTdphi(:,:,k);
                    dSigma_timedphi(:,k) = temp(ind_t);
                end
            end
            
            logL_D = logL_D - J_D;
            logL_T = logL_T - J_T;
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
                
                dJ_Ddphi = chainrule_dxdy_dydz(dJ_DdY,dY_sidphi) + chainrule_dxdy_dydz(dJ_DdSigma,dSigma_noisedphi) ;
                dJ_Ddxi = chainrule_dxdy_dydz(dJ_Ddphi,dphidxi);
                
                dJ_Tdphi = chainrule_dxdy_dydz(dJ_TdT,dT_sidphi) + chainrule_dxdy_dydz(dJ_TdSigma,dSigma_timedphi) ;
                dJ_Tdxi = chainrule_dxdy_dydz(dJ_Tdphi,dphidxi);
                
                dbdxi = chainrule_dxdy_dydz(dbhat_sidbeta,dbetadxi) + chainrule_dxdy_dydz(dbhat_siddelta,ddeltadxi);
                P{s}.SCTL.dbdxi(:,:,i) = dbdxi;
                
                dJ_bdxi = chainrule_dxdy_dydz(dJ_bdb,dbdxi) + chainrule_dxdy_dydz(pdJ_bpddelta,ddeltadxi);
                
                dlogLdxi = dlogLdxi - transpose(dJ_Ddxi) - transpose(dJ_Tdxi) - transpose(dJ_bdxi);
                
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
                        bsxfun(@times,ddJ_DdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
                    
                    ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dY_sidphi,[3,1,2])) + ...
                        bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
                    
                    ddJ_Ddphidphi = chainrule_dxdy_dydz(dJ_DdY,ddY_sidphidphi) ...
                        + squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(dY_sidphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddJ_DdphidSigma,permute(dSigma_noisedphi,[3,1,4,2])),2));
                    
                    ddJ_Ddxidxi = chainrule_dxdy_dydz(dJ_Ddphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Ddphidphi,dphidxi);
                    
                    ddJ_TdphidT = bsxfun(@times,ddJ_TdTdT,permute(dT_sidphi,[3,1,2])) + ...
                        bsxfun(@times,ddJ_TdTdSigma,permute(dSigma_timedphi,[3,1,2]));
                    
                    ddJ_TdphidSigma = bsxfun(@times,ddJ_TdTdSigma,permute(dT_sidphi,[3,1,2])) + ...
                        bsxfun(@times,ddJ_TdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    
                    ddJ_Tdphidphi = chainrule_dxdy_dydz(dJ_TdT,ddT_sidphidphi) ...
                        + squeeze(sum(bsxfun(@times,ddJ_TdphidT,permute(dT_sidphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddJ_TdphidSigma,permute(dSigma_timedphi,[3,1,4,2])),2));
                    
                    ddJ_Tdxidxi = chainrule_dxdy_dydz(dJ_Tdphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Tdphidphi,dphidxi);
                    
                    ddbdxidxi = chainrule_dxdy_dydz(dbhat_sidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddbhat_sidbetadbeta,dbetadxi) ...
                        + chainrule_dxdy_dydz(dbhat_siddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddbhat_siddeltaddelta,ddeltadxi) ...
                        + chainrule_ddxdydy_dydz_dydv(ddbhat_sidbetaddelta,dbetadxi,ddeltadxi) ...
                        + chainrule_ddxdydy_dydz_dydv(permute(ddbhat_sidbetaddelta,[1,3,2]),ddeltadxi,dbetadxi);
                    
                    P{s}.SCTL.ddbdxidxi(:,:,:,i) = ddbdxidxi;
                    
                    ddJ_bdxidxi = chainrule_dxdy_dydz(dJ_bdb,ddbdxidxi) + chainrule_ddxdydy_dydz(ddJ_bdbdb,dbdxi) ...
                        + chainrule_dxdy_dydz(pdJ_bpddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(pdpdJ_bpddeltapddelta,ddeltadxi) ...
                        + chainrule_ddxdydy_dydz_dydv(dpdJ_bdbpddelta,dbdxi,ddeltadxi) ...
                        + chainrule_ddxdydy_dydz_dydv(permute(dpdJ_bdbpddelta,[2,1]),ddeltadxi,dbdxi);
                    
                    ddlogLdxidxi = ddlogLdxidxi - ddJ_Ddxidxi - ddJ_Tdxidxi - ddJ_bdxidxi;
                    
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
            Sim_SCTL.Y(:,:,i) = Y(ind_time,:);
            Sim_SCTL.T(:,:,i) = T(ind_t,:);
            Sim_SCTL.R(:,:,i) = R(ind_t,:);
        end
        
        logL = logL + logL_D + logL_T + logL_b + logL_I;
        
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
                    bar([logL_D,logL_T,logL_b,logL_I,logL_s])
                    set(gca,'XTickLabel',{'Data','Event','Par','Int','Pen'})
                else
                    bar([logL_D,logL_T,logL_b,logL_s])
                    set(gca,'XTickLabel',{'Data','Event','Par','Pen'})
                end
            else
                if(Model.integration)
                    bar([logL_D,logL_T,logL_b,logL_I])
                    set(gca,'XTickLabel',{'Data','Event','Par','Int'})
                else
                    bar([logL_D,logL_T,logL_b])
                    set(gca,'XTickLabel',{'Data','Event','Par'})
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
function varargout = optimize_SCTL_si(Model,Data,bhat_0,beta,delta,type_D,t,Ym,Tm,ind_y,ind_t,F_diff,b_diff,s)
options_fmincon = optimset('algorithm','trust-region-reflective',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',1000,... % 1000%'display','iter',...
    'TolFun',1e-6,...
    'TolX',1e-6,...
    'Hessian','user-supplied');

[bhat,~,~,~,~,~,~] = fmincon(...
    @(b) objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),...
    bhat_0,[],[],[],[],-10*ones(length(bhat_0),1),10*ones(length(bhat_0),1),[],options_fmincon);


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
                [~,~,G,...
                    ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta,...
                    dGdb,pdGpdbeta,pdGpddelta] = ...
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
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
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-4,2,3)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(beta,@(beta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-4,2,4)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(delta)objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-4,2,5)
                % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat,@(b)objective_SCTL_s1(Model,beta,b,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s),1e-4,1,2)
                [~,~,G,...
                    ddJdbdbeta,ddJdbddelta,ddJdbetadbeta,ddJddeltaddelta,ddJdbetaddelta] = ...
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
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
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind,s);
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
                    objective_SCTL_s1(Model,beta,bhat,Data{s}.condition,delta,type_D,t,Ym,Tm,ind_y,ind_t,s);
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

function varargout = objective_SCTL_s1(Model,beta,b,kappa,delta,type_D,t,Ym,Tm,ind_y,ind_t,s)

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);

% Single-cell parameters
phi = Model.exp{s}.phi(beta,b);

sigma_noise = Model.exp{s}.sigma_noise(phi);

if(size(sigma_noise,1) == size(Ym,1))
    if(size(sigma_noise,2) == 1)
        Sigma_noise = repmat(sigma_noise,[1,size(Ym,2)]);
    elseif(size(sigma,2) == size(Ym,2))
        Sigma_noise = sigma_noise;
    else
        error('Incompatible size of sigma_noise parametrisation!')
    end
elseif(size(sigma_noise,2) == size(Ym,2))
    if(size(sigma_noise,1) == 1)
        Sigma_noise = repmat(sigma_noise,[size(Ym,1),1]);
    else
        error('Incompatible size of sigma_noise parametrisation!')
    end
elseif(and(size(sigma_noise,1)==1,size(sigma_noise,2)==1))
    Sigma_noise = repmat(sigma_noise,size(Ym));
else
    error('Incompatible size of sigma_noise parametrisation!')
end

sigma_time = Model.exp{s}.sigma_time(phi);

if(size(sigma_time,1) == size(Tm,1))
    if(size(sigma_time,2) == 1)
        Sigma_time = repmat(sigma_time,[1,size(Tm,2)]);
    elseif(size(sigma,2) == size(Tm,2))
        Sigma_time = sigma_time;
    else
        error('Incompatible size of sigma_time parametrisation!')
    end
elseif(size(sigma_time,2) == size(Tm,2))
    if(size(sigma_time,1) == 1)
        Sigma_time = repmat(sigma_time,[size(Tm,1),1]);
    else
        error('Incompatible size of sigma_time parametrisation!')
    end
elseif(and(size(sigma_time,1)==1,size(sigma_time,2)==1))
    Sigma_time = repmat(sigma_time,size(Tm));
else
    error('Incompatible size of sigma_time parametrisation!')
end



% Simulate model
if nargout >= 3
    %[~,~,~,Y,~,dYdphi,~,ddYdphidphi] = Model.exp{s}.model(t,phi,kappa);
    option_simu.sensi = 1;
    sol = Model.exp{s}.model(t,phi,kappa,option_simu);
    Y = sol.y;
    dYdphi = sol.sy;
    T = sol.root(ind_t,:);
    R = sol.rootval(ind_t,:);
    temp1 = dYdphi;
    temp2 = sol.rootS(ind_t,:,:);
    %         temp2 = ddYdphidphi;
    dYdphi = zeros(length(ind_y),length(phi));
    dTdphi = zeros(length(ind_t),length(phi));
    ddYdphidphi = zeros(length(ind_y),length(phi),length(phi));
    ddTdphidphi = zeros(length(ind_t),length(phi),length(phi));
    for j = 1:length(phi)
        tempy = temp1(:,:,j);
        tempt = temp2(:,:,j);
        dYdphi(:,j) = tempy(ind_y);
        dTdphi(:,j) = tempt(ind_t);
        %             for k = 1:length(phi)
        %                 temp = temp2(:,:,j,k);
        %                 ddYdphidphi(:,j,k) = temp(ind);
        %             end
    end
elseif nargout >=2
    option_simu.sensi = 1;
    sol = Model.exp{s}.model(t,phi,kappa,option_simu);
    Y = sol.y;
    dYdphi = sol.sy;
    T = sol.root(ind_t,:);
    R = sol.rootval(ind_t,:);
    temp1 = dYdphi;
    temp2 = sol.rootS(ind_t,:,:);
    dYdphi = zeros(length(ind_y),length(phi));
    dTdphi = zeros(length(ind_t),length(phi));
    for j = 1:length(phi)
        tempy = temp1(:,:,j);
        tempt = temp2(:,:,j);
        dYdphi(:,j) = tempy(ind_y);
        dTdphi(:,j) = tempt(ind_t);
    end
else
    option_simu.sensi = 0;
    sol = Model.exp{s}.model(t,phi,kappa,option_simu);
    Y = sol.y;
    T = sol.root(ind_t,:);
    R = sol.rootval(ind_t,:);
end

% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)Model.exp{s}.model(t,phi,kappa),1e-4,4,6)

if nargout >= 2 % first order derivatives
    dsigma_noisedphi = Model.exp{s}.dsigma_noisedphi(phi);
    dSigma_noisedphi = zeros(length(ind_y),length(phi));
    
    if(size(dsigma_noisedphi,1) == size(Ym,1))
        if(size(dsigma_noisedphi,2) == 1)
            dSNdphi = repmat(dsigma_noisedphi,[1,size(Ym,2),1]);
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            dSNdphi = dsigma_noisedphi;
        end
    elseif(size(dsigma_noisedphi,2) == size(Ym,2))
        if(size(dsigma_noisedphi,1) == 1)
            dSNdphi = repmat(dsigma_noisedphi,[size(Ym,1),1,1]);
        end
    elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
        dSNdphi = repmat(dsigma_noisedphi,[size(Ym),1]);
    end
    
    if nargout >= 3 % second order derivatives
        ddsigma_noisedphidphi = Model.exp{s}.ddsigma_noisedphidphi(phi);
        ddSigma_noisedphidphi = zeros(length(ind_y),length(phi),length(phi));
        if(size(dsigma_noisedphi,1) == size(Ym,1))
            if(size(dsigma_noisedphi,2) == 1)
                ddSNdphidphi = repmat(ddsigma_noisedphidphi,[1,size(Ym,2),1,1]);
            elseif(size(dsigma_noisedphi,2) == size(Ym,2))
                ddSNdphidphi = ddsigma_noisedphidphi ;
            end
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            if(size(dsigma_noisedphi,1) == 1)
                ddSNdphidphi = repmat(ddsigma_noisedphidphi,[size(Ym,1),1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            ddSNdphidphi = repmat(ddsigma_noisedphidphi,[size(Ym),1,1]);
        end
    end
    
    if nargout >= 9 % third order derivatives
        dddsigma_noisedphidphidphi = Model.exp{s}.dddsigma_noisedphidphidphi(phi);
        dddSigma_noisedphidphidphi = zeros(length(ind_y),length(phi),length(phi),length(phi));
        if(size(dsigma_noisedphi,1) == size(Ym,1))
            if(size(dsigma_noisedphi,2) == 1)
                dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[1,size(Ym,2),1,1,1]);
            elseif(size(dsigma_noisedphi,2) == size(Ym,2))
                dddSNdphidphidphi = dddsigma_noisedphidphidphi;
            end
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            if(size(dsigma_noisedphi,1) == 1)
                dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[size(Ym,1),1,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[size(Ym),1,1,1]);
        end
    end
    
    if nargout >= 14 % fourth order derivatives
        ddddsigma_noisedphidphidphidphi = Model.exp{s}.ddddsigma_noisedphidphidphidphi(phi);
        ddddSigma_noisedphidphidphidphi = zeros(length(ind_y),length(phi),length(phi),length(phi),length(phi));
        if(size(dsigma_noisedphi,1) == size(Ym,1))
            if(size(dsigma_noisedphi,2) == 1)
                ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[1,size(Ym,2),1,1,1]);
            elseif(size(dsigma_noisedphi,2) == size(Ym,2))
                ddddSNdphidphidphidphi = ddddsigma_noisedphidphidphidphi;
            end
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            if(size(dsigma_noisedphi,1) == 1)
                ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[size(Ym,1),1,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[size(Ym),1,1,1]);
        end
    end
    
    for k = 1:length(phi) % first order derivatives
        temp = dSNdphi(:,:,k);
        dSigma_noisedphi(:,k) = temp(ind_y);
        if nargout >= 3 % second order derivatives
            for l = 1:length(phi)
                temp = ddSNdphidphi(:,:,k,l);
                ddSigma_noisedphidphi(:,k,l) = temp(ind_y);
                if nargout >= 9 % third order derivatives
                    for m = 1:length(phi)
                        temp = dddSNdphidphidphi(:,:,k,l,m);
                        dddSigma_noisedphidphidphi(:,k,l,m) = temp(ind_y);
                        if nargout >= 14 % fourth order derivatives
                            for n = 1:length(phi)
                                temp = ddddSNdphidphidphidphi(:,:,k,l,m,n);
                                ddddSigma_noisedphidphidphidphi(:,k,l,m,n) = temp(ind_y);
                            end
                        end
                    end
                end
            end
        end
    end
end

if nargout >= 2 % first order derivatives
    dsigma_timedphi = Model.exp{s}.dsigma_timedphi(phi);
    dSigma_timedphi = zeros(length(ind_t),length(phi));
    
    if(size(dsigma_timedphi,1) == size(Tm,1))
        if(size(dsigma_timedphi,2) == 1)
            dSTdphi = repmat(dsigma_timedphi,[1,size(Tm,2),1]);
        elseif(size(dsigma_timedphi,2) == size(Tm,2))
            dSTdphi = dsigma_timedphi;
        end
    elseif(size(dsigma_timedphi,2) == size(Tm,2))
        if(size(dsigma_timedphi,1) == 1)
            dSTdphi = repmat(dsigma_timedphi,[size(Tm,1),1,1]);
        end
    elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
        dSTdphi = repmat(dsigma_timedphi,[size(Tm),1]);
    end
    
    if nargout >= 3 % second order derivatives
        ddsigma_timedphidphi = Model.exp{s}.ddsigma_timedphidphi(phi);
        ddSigma_timedphidphi = zeros(length(ind_t),length(phi),length(phi));
        if(size(dsigma_timedphi,1) == size(Tm,1))
            if(size(dsigma_timedphi,2) == 1)
                ddSTdphidphi = repmat(ddsigma_timedphidphi,[1,size(Tm,2),1,1]);
            elseif(size(dsigma_timedphi,2) == size(Tm,2))
                ddSTdphidphi = ddsigma_timedphidphi ;
            end
        elseif(size(dsigma_timedphi,2) == size(Tm,2))
            if(size(dsigma_timedphi,1) == 1)
                ddSTdphidphi = repmat(ddsigma_timedphidphi,[size(Tm,1),1,1]);
            end
        elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
            ddSTdphidphi = repmat(ddsigma_timedphidphi,[size(Tm),1,1]);
        end
    end
    
    if nargout >= 9 % third order derivatives
        dddsigma_timedphidphidphi = Model.exp{s}.dddsigma_timedphidphidphi(phi);
        dddSigma_timedphidphidphi = zeros(length(ind_t),length(phi),length(phi),length(phi));
        if(size(dsigma_timedphi,1) == size(Tm,1))
            if(size(dsigma_timedphi,2) == 1)
                dddSTdphidphidphi = repmat(dddsigma_timedphidphidphi,[1,size(Tm,2),1,1,1]);
            elseif(size(dsigma_timedphi,2) == size(Tm,2))
                dddSTdphidphidphi = dddsigma_timedphidphidphi;
            end
        elseif(size(dsigma_timedphi,2) == size(Tm,2))
            if(size(dsigma_timedphi,1) == 1)
                dddSTdphidphidphi = repmat(dddsigma_timedphidphidphi,[size(Tm,1),1,1,1]);
            end
        elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
            dddSTdphidphidphi = repmat(dddsigma_timedphidphidphi,[size(Tm),1,1,1]);
        end
    end
    
    if nargout >= 14 % fourth order derivatives
        ddddsigma_timedphidphidphidphi = Model.exp{s}.ddddsigma_timedphidphidphidphi(phi);
        ddddSigma_timedphidphidphidphi = zeros(length(ind_t),length(phi),length(phi),length(phi),length(phi));
        if(size(dsigma_timedphi,1) == size(Tm,1))
            if(size(dsigma_timedphi,2) == 1)
                ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[1,size(Tm,2),1,1,1]);
            elseif(size(dsigma_timedphi,2) == size(Tm,2))
                ddddSTdphidphidphidphi = ddddsigma_timedphidphidphidphi;
            end
        elseif(size(dsigma_timedphi,2) == size(Tm,2))
            if(size(dsigma_timedphi,1) == 1)
                ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[size(Tm,1),1,1,1]);
            end
        elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
            ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[size(Tm),1,1,1]);
        end
    end
    
    for k = 1:length(phi) % first order derivatives
        temp = dSTdphi(:,:,k);
        dSigma_timedphi(:,k) = temp(ind_t);
        if nargout >= 3 % second order derivatives
            for l = 1:length(phi)
                temp = ddSTdphidphi(:,:,k,l);
                ddSigma_timedphidphi(:,k,l) = temp(ind_t);
                if nargout >= 9 % third order derivatives
                    for m = 1:length(phi)
                        temp = dddSTdphidphidphi(:,:,k,l,m);
                        dddSigma_timedphidphidphi(:,k,l,m) = temp(ind_t);
                        if nargout >= 14 % fourth order derivatives
                            for n = 1:length(phi)
                                temp = ddddSTdphidphidphidphi(:,:,k,l,m,n);
                                ddddSigma_timedphidphidphidphi(:,k,l,m,n) = temp(ind_t);
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
            J_D = normal_noise(Y,Ym,Sigma_noise,ind_y);
        elseif nargout <=2 % first order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma] = normal_noise(Y,Ym,Sigma_noise,ind_y);
        elseif nargout <=8 % second order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = normal_noise(Y,Ym,Sigma_noise,ind_y);
        elseif nargout <=13 % third order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma] = normal_noise(Y,Ym,Sigma_noise,ind_y);
        else % fourth order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma,...
                ddddJ_DdYdYdYdY,ddddJ_DdYdYdYdSigma,ddddJ_DdYdYdSigmadSigma,ddddJ_DdYdSigmadSigmadSigma,ddddJ_DdSigmadSigmadSigmadSigma] = normal_noise(Y,Ym,Sigma_noise,ind_y);
        end
        
    case 'lognormal'
        
        if nargout <=1
            J_D = lognormal_noise(Y,Ym,Sigma_noise,ind_y);
        elseif nargout <=2 % first order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma] = lognormal_noise(Y,Ym,Sigma_noise,ind_y);
        elseif nargout <=8 % second order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = lognormal_noise(Y,Ym,Sigma_noise,ind_y);
        elseif nargout <=13 % third order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma] = lognormal_noise(Y,Ym,Sigma_noise,ind_y);
        else % fourth order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma,...
                ddddJ_DdYdYdYdY,ddddJ_DdYdYdYdSigma,ddddJ_DdYdYdSigmadSigma,ddddJ_DdYdSigmadSigmadSigma,ddddJ_DdSigmadSigmadSigmadSigma] = lognormal_noise(Y,Ym,Sigma_noise,ind_y);
        end
        
end

% time model
switch(Model.exp{s}.time_model)
    case 'normal'
        if nargout <= 1
            J_T = normal_time(T,Tm,R,Sigma_time,ind_t);
        elseif nargout <=2 % first order derivatives
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(T,@(T)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,1,2)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,1,3)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(T,@(T)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,2,4)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,2,5)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,3,6)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(T,@(T)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,4,7)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,4,8)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,5,9)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,6,10)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(T,@(T)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,7,11)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,7,12)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,8,13)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,9,14)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,Tm,R,Sigma_time,ind_t),1e-5,10,15)
            [J_T,...
                dJ_TdT,dJ_TdSigma] = normal_time(T,Tm,R,Sigma_time,ind_t);
        elseif nargout <=8 % second order derivatives
            [J_T,...
                dJ_TdT,dJ_TdSigma,...
                ddJ_TdTdT,ddJ_TdTdSigma,ddJ_TdSigmadSigma] = normal_time(T,Tm,R,Sigma_time,ind_t);
        elseif nargout <=13 % third order derivatives
            [J_T,...
                dJ_TdT,dJ_TdSigma,...
                ddJ_TdTdT,ddJ_TdTdSigma,ddJ_TdSigmadSigma,...
                dddJ_TdTdTdT,dddJ_TdTdTdSigma,dddJ_TdTdSigmadSigma,dddJ_TdSigmadSigmadSigma] = normal_time(T,Tm,R,Sigma_time,ind_t);
        else % fourth order derivatives
            [J_T,...
                dJ_TdT,dJ_TdSigma,...
                ddJ_TdTdT,ddJ_TdTdSigma,ddJ_TdSigmadSigma,...
                dddJ_TdTdTdT,dddJ_TdTdTdSigma,dddJ_TdTdSigmadSigma,dddJ_TdSigmadSigmadSigma,...
                ddddJ_TdTdTdTdT,ddddJ_TdTdTdTdSigma,ddddJ_TdTdTdSigmadSigma,ddddJ_TdTdSigmadSigmadSigma,ddddJ_TdSigmadSigmadSigmadSigma] = normal_time(T,Tm,R,Sigma_time,ind_t);
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

J = J_D + J_T + J_b ;
varargout{1} = J;

if nargout >= 2
    %% dJdb
    dphidb = Model.exp{s}.dphidb(beta,b);
    
    dJ_Ddphi = chainrule_dxdy_dydz(dJ_DdY,dYdphi) + chainrule_dxdy_dydz(dJ_DdSigma,dSigma_noisedphi);
    
    dJ_Tdphi = chainrule_dxdy_dydz(dJ_TdT,dTdphi) + chainrule_dxdy_dydz(dJ_TdSigma,dSigma_timedphi);
    
    dJ_Ddb = chainrule_dxdy_dydz(dJ_Ddphi,dphidb);
    
    dJ_Tdb = chainrule_dxdy_dydz(dJ_Tdphi,dphidb);
    
    dJdb = dJ_Ddb + dJ_Tdb + dJ_bdb;
    
    varargout{2} = dJdb;
    
    if nargout >= 3
        %% ddJdbdb
        
        ddJ_DdphidY = bsxfun(@times,ddJ_DdYdY,permute(dYdphi,[3,1,2])) + bsxfun(@times,ddJ_DdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
        ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dYdphi,[3,1,2])) + bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
        ddJ_Ddphidphi = squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(dYdphi,[3,1,4,2])) + bsxfun(@times,ddJ_DdphidSigma,permute(dSigma_noisedphi,[3,1,4,2])),2));
        
        ddphidbdb = Model.exp{s}.ddphidbdb(beta,b);
        
        ddJ_Ddbdphi = transpose(squeeze(sum(bsxfun(@times,ddJ_Ddphidphi,permute(dphidb,[3,1,2])),2)));
        ddJ_Ddbdb = squeeze(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidb,[3,1,2,4])),2)) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbdb);
        
        ddJ_TdphidT = bsxfun(@times,ddJ_TdTdT,permute(dTdphi,[3,1,2])) + bsxfun(@times,ddJ_TdTdSigma,permute(dSigma_timedphi,[3,1,2]));
        ddJ_TdphidSigma = bsxfun(@times,ddJ_TdTdSigma,permute(dTdphi,[3,1,2])) + bsxfun(@times,ddJ_TdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
        ddJ_Tdphidphi = squeeze(sum(bsxfun(@times,ddJ_TdphidT,permute(dTdphi,[3,1,4,2])) + bsxfun(@times,ddJ_TdphidSigma,permute(dSigma_timedphi,[3,1,4,2])),2));
        
        ddphidbdb = Model.exp{s}.ddphidbdb(beta,b);
        
        ddJ_Tdbdphi = transpose(squeeze(sum(bsxfun(@times,ddJ_Tdphidphi,permute(dphidb,[3,1,2])),2)));
        ddJ_Tdbdb = squeeze(sum(bsxfun(@times,ddJ_Tdbdphi,permute(dphidb,[3,1,2,4])),2)) + chainrule_dxdy_dydz(dJ_Tdphi,ddphidbdb);
        
        
        ddJdbdb = squeeze(ddJ_Ddbdb) + squeeze(ddJ_Tdbdb) + squeeze(ddJ_bdbdb);
        
        varargout{3} = ddJdbdb;
        
        if nargout >= 4
            
            %% ddJdbdbeta
            dphidbeta = Model.exp{s}.dphidbeta(beta,b);
            ddphidbdbeta = Model.exp{s}.ddphidbdbeta(beta,b);
            
            % if size of b == 1 then we have to permute second term,
            if(numel(b) == 1)
                ddJ_Ddbdbeta = permute(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidbeta,[3,1,2])),2),[1,3,4,2]) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbdbeta);
                ddJ_Tdbdbeta = permute(sum(bsxfun(@times,ddJ_Tdbdphi,permute(dphidbeta,[3,1,2])),2),[1,3,4,2]) + chainrule_dxdy_dydz(dJ_Tdphi,ddphidbdbeta);
            else
                ddJ_Ddbdbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidbeta,[3,1,2])),2)) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbdbeta);
                ddJ_Tdbdbeta = squeeze(sum(bsxfun(@times,ddJ_Tdbdphi,permute(dphidbeta,[3,1,2])),2)) + chainrule_dxdy_dydz(dJ_Tdphi,ddphidbdbeta);
            end
            
            ddJdbdbeta = ddJ_Ddbdbeta + ddJ_Tdbdbeta;
            
            varargout{4} = ddJdbdbeta;
            
            %% ddJdbdelta
            
            ddJdbddelta = ddJ_bdbddelta;
            
            varargout{5} = ddJdbddelta;
            
            %% ddJdbetadbeta
            
            ddphidbetadbeta = Model.exp{s}.ddphidbetadbeta(beta,b);
            
            ddJ_Ddbetadphi = permute(sum(bsxfun(@times,ddJ_Ddphidphi,permute(dphidbeta,[3,1,2])),2),[3,2,1]);
            ddJ_Tdbetadphi = permute(sum(bsxfun(@times,ddJ_Tdphidphi,permute(dphidbeta,[3,1,2])),2),[3,2,1]);
            if(numel(b)==1)
                
            else
                ddJ_Ddbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule_dxdy_dydz(dJ_Ddphi,ddphidbetadbeta);
                ddJ_Tdbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Tdbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule_dxdy_dydz(dJ_Tdphi,ddphidbetadbeta);
            end
            
            ddJdbetadbeta = ddJ_Ddbetadbeta + ddJ_Tdbetadbeta;
            
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
                    + bsxfun(@times,ddJ_DdphidSigma,permute(ddSigma_noisedphidphi,[4,1,5,2,3])),2));
                
                dddJ_Ddphidphidphi = chainrule_dxdy_dydz(dJ_DdSigma,dddSigma_noisedphidphidphi) ...
                    + temp + permute(temp,[2,1,3]) + permute(temp,[2,3,1]);
                
                dddJ_DdphidYdY = bsxfun(@times,dddJ_DdYdYdY,permute(dYdphi,[3,1,2])) ...
                    + bsxfun(@times,dddJ_DdYdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
                dddJ_DdphidYdSigma = bsxfun(@times,dddJ_DdYdYdSigma,permute(dYdphi,[3,1,2])) ...
                    + bsxfun(@times,dddJ_DdYdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
                dddJ_DdphidSigmadSigma = bsxfun(@times,dddJ_DdYdSigmadSigma,permute(dYdphi,[3,1,2])) ...
                    + bsxfun(@times,dddJ_DdSigmadSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
                dddJ_DdphidphidY = bsxfun(@times,dddJ_DdphidYdY,permute(dYdphi,[3,1,4,2])) ...
                    + bsxfun(@times,dddJ_DdphidYdSigma,permute(dSigma_noisedphi,[3,1,4,2]));
                dddJ_DdphidphidSigma = bsxfun(@times,dddJ_DdphidYdSigma,permute(dYdphi,[3,1,4,2])) ...
                    + bsxfun(@times,dddJ_DdphidSigmadSigma,permute(dSigma_noisedphi,[3,1,4,2]));
                
                dddJ_Ddphidphidphi = dddJ_Ddphidphidphi + squeeze(sum(bsxfun(@times,dddJ_DdphidphidY,permute(dYdphi,[3,1,4,5,2])) ...
                    + bsxfun(@times,dddJ_DdphidphidSigma,permute(dSigma_noisedphi,[3,1,4,5,2])),2));
                
                dddJ_Ddbdphidphi = permute(sum(bsxfun(@times,dddJ_Ddphidphidphi,permute(dphidb,[1,3,4,2])),1),[4,2,3,1]);
                dddJ_Ddbdbdphi = permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(dphidb,[3,1,4,2])),2),[1,4,3,2]) ...
                    + permute(sum(bsxfun(@times,ddJ_Ddphidphi,permute(ddphidbdb,[1,4,2,3])),1),[3,4,2,1]);
                dddJ_Ddbdbdb = permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(dphidb,[3,4,1,2])),3),[1,2,4,3]) ...
                    + squeeze(chainrule_dxdy_dydz(ddJ_Ddbdphi,ddphidbdb));
                
                temp = squeeze(sum(bsxfun(@times,ddJ_TdphidT,permute(ddTdphidphi,[4,1,5,2,3])) ...
                    + bsxfun(@times,ddJ_TdphidSigma,permute(ddSigma_timedphidphi,[4,1,5,2,3])),2));
                
                dddJ_Tdphidphidphi = chainrule_dxdy_dydz(dJ_TdSigma,dddSigma_timedphidphidphi) ...
                    + temp + permute(temp,[2,1,3]) + permute(temp,[2,3,1]);
                
                
                dddJ_TdphidTdT = bsxfun(@times,dddJ_TdTdTdT,permute(dTdphi,[3,1,2])) ...
                    + bsxfun(@times,dddJ_TdTdTdSigma,permute(dSigma_timedphi,[3,1,2]));
                dddJ_TdphidTdSigma = bsxfun(@times,dddJ_TdTdTdSigma,permute(dTdphi,[3,1,2])) ...
                    + bsxfun(@times,dddJ_TdTdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                dddJ_TdphidSigmadSigma = bsxfun(@times,dddJ_TdTdSigmadSigma,permute(dTdphi,[3,1,2])) ...
                    + bsxfun(@times,dddJ_TdSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                dddJ_TdphidphidT = bsxfun(@times,dddJ_TdphidTdT,permute(dTdphi,[3,1,4,2])) ...
                    + bsxfun(@times,dddJ_TdphidTdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                dddJ_TdphidphidSigma = bsxfun(@times,dddJ_TdphidTdSigma,permute(dTdphi,[3,1,4,2])) ...
                    + bsxfun(@times,dddJ_TdphidSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                
                dddJ_Tdphidphidphi = dddJ_Tdphidphidphi + squeeze(sum(bsxfun(@times,dddJ_TdphidphidT,permute(dTdphi,[3,1,4,5,2])) ...
                    + bsxfun(@times,dddJ_TdphidphidSigma,permute(dSigma_timedphi,[3,1,4,5,2])),2));
                
                dddJ_Tdbdphidphi = permute(sum(bsxfun(@times,dddJ_Tdphidphidphi,permute(dphidb,[1,3,4,2])),1),[4,2,3,1]);
                dddJ_Tdbdbdphi = permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(dphidb,[3,1,4,2])),2),[1,4,3,2]) ...
                    + permute(sum(bsxfun(@times,ddJ_Tdphidphi,permute(ddphidbdb,[1,4,2,3])),1),[3,4,2,1]);
                dddJ_Tdbdbdb = permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(dphidb,[3,4,1,2])),3),[1,2,4,3]) ...
                    + squeeze(chainrule_dxdy_dydz(ddJ_Tdbdphi,ddphidbdb));
                
                dddJdbdbdb = dddJ_Ddbdbdb + dddJ_Tdbdbdb + squeeze(dddJ_bdbdbdb);
                
                varargout{9} = dddJdbdbdb;
                
                %% dddJdbdbbeta
                
                dddJ_Ddbdbdbeta = permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                dddJ_Tdbdbdbeta = permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                dddJdbdbdbeta = dddJ_Ddbdbdbeta + dddJ_Tdbdbdbeta;
                
                varargout{10} = dddJdbdbdbeta;
                
                
                %% dddJdbdbdelta
                dddJdbdbddelta = dddJ_bdbdbddelta;
                
                varargout{11} = squeeze(dddJdbdbddelta);
                
                %% dddJdbdbetadbeta
                
                dddJ_Ddbdbetadphi = permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(dphidbeta,[3,1,4,2])),2),[1,4,3,2]);
                dddJ_Ddbdbetadbeta = permute(sum(bsxfun(@times,dddJ_Ddbdbetadphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                dddJ_Tdbdbetadphi = permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(dphidbeta,[3,1,4,2])),2),[1,4,3,2]);
                dddJ_Tdbdbetadbeta = permute(sum(bsxfun(@times,dddJ_Tdbdbetadphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                dddJdbdbetadbeta = dddJ_Ddbdbetadbeta + dddJ_Tdbdbetadbeta;
                
                varargout{12} = squeeze(dddJdbdbetadbeta);
                
                %% dddJdbdbdelta
                dddJdbddeltaddelta = dddJ_bdbddeltaddelta;
                
                varargout{13} = squeeze(dddJdbddeltaddelta);
                
                
                if nargout >= 14
                    
                    %% ddddJdbdbdbdb
                    
                    temp = squeeze(sum(bsxfun(@times,bsxfun(@times,ddJ_DdYdY,permute(ddYdphidphi,[4,1,2,3])),permute(ddYdphidphi,[4,1,5,6,2,3])) ...
                        + bsxfun(@times,bsxfun(@times,ddJ_DdYdSigma,permute(ddYdphidphi,[4,1,2,3])),permute(ddSigma_noisedphidphi,[4,1,5,6,2,3])) ...
                        + bsxfun(@times,bsxfun(@times,ddJ_DdSigmadSigma,permute(ddSigma_noisedphidphi,[4,1,2,3])),permute(ddSigma_noisedphidphi,[4,1,5,6,2,3])),2));
                    
                    ddddJ_Ddphidphidphidphi = ...
                        permute(temp,[1,3,2,4]) + permute(temp,[1,3,4,2]) + permute(temp,[3,4,1,2]);
                    
                    ddddJ_DdphidYdYdY = bsxfun(@times,ddddJ_DdYdYdYdY,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_DdYdYdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
                    ddddJ_DdphidYdYdSigma = bsxfun(@times,ddddJ_DdYdYdYdSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_DdYdYdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
                    ddddJ_DdphidYdSigmadSigma = bsxfun(@times,ddddJ_DdYdYdSigmadSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_DdYdSigmadSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
                    ddddJ_DdphidSigmadSigmadSigma = bsxfun(@times,ddddJ_DdYdSigmadSigmadSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_DdSigmadSigmadSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
                    ddddJ_DdphidphidYdY = bsxfun(@times,ddddJ_DdphidYdYdY,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddddJ_DdphidYdYdSigma,permute(dSigma_noisedphi,[3,1,4,2]));
                    ddddJ_DdphidphidYdSigma = bsxfun(@times,ddddJ_DdphidYdYdSigma,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddddJ_DdphidYdSigmadSigma,permute(dSigma_noisedphi,[3,1,4,2]));
                    ddddJ_DdphidphidSigmadSigma = bsxfun(@times,ddddJ_DdphidYdSigmadSigma,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddddJ_DdphidSigmadSigmadSigma,permute(dSigma_noisedphi,[3,1,4,2]));
                    ddddJ_DdphidphidphidY = bsxfun(@times,ddddJ_DdphidphidYdY,permute(dYdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,ddddJ_DdphidphidYdSigma,permute(dSigma_noisedphi,[3,1,4,5,2]));
                    ddddJ_DdphidphidphidSigma = bsxfun(@times,ddddJ_DdphidphidYdSigma,permute(dYdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,ddddJ_DdphidphidSigmadSigma,permute(dSigma_noisedphi,[3,1,4,5,2]));
                    ddddJ_Ddphidphidphidphi = ddddJ_Ddphidphidphidphi + squeeze(sum(bsxfun(@times,ddddJ_DdphidphidphidY,permute(dYdphi,[3,1,4,5,6,2])) ...
                        + bsxfun(@times,ddddJ_DdphidphidphidSigma,permute(dSigma_noisedphi,[3,1,4,5,6,2])),2));
                    
                    ddddJ_Ddbdphidphidphi = permute(sum(bsxfun(@times,ddddJ_Ddphidphidphidphi,permute(dphidb,[1,3,4,5,2])),1),[5,2,3,4,1]);
                    ddddJ_Ddbdbdphidphi = permute(sum(bsxfun(@times,ddddJ_Ddbdphidphidphi,permute(dphidb,[3,1,4,5,2])),2),[1,5,3,4,2]) ...
                        + permute(sum(bsxfun(@times,dddJ_Ddphidphidphi,permute(ddphidbdb,[1,4,5,2,3])),1),[4,5,2,3,1]);
                    ddddJ_Ddbdbdbdphi = permute(sum(bsxfun(@times,ddddJ_Ddbdbdphidphi,permute(dphidb,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                        + permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(ddphidbdb,[4,1,5,2,3])),2),[1,4,5,3,2]);
                    ddddJ_Ddbdbdbdb = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbdphi,permute(dphidb,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                        + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbdb,[4,5,1,2,3])),3),[1,2,4,5,3]);
                    
                    temp = squeeze(sum(bsxfun(@times,bsxfun(@times,ddJ_TdTdT,permute(ddYdphidphi,[4,1,2,3])),permute(ddYdphidphi,[4,1,5,6,2,3])) ...
                        + bsxfun(@times,bsxfun(@times,ddJ_TdTdSigma,permute(ddYdphidphi,[4,1,2,3])),permute(ddSigma_timedphidphi,[4,1,5,6,2,3])) ...
                        + bsxfun(@times,bsxfun(@times,ddJ_TdSigmadSigma,permute(ddSigma_timedphidphi,[4,1,2,3])),permute(ddSigma_timedphidphi,[4,1,5,6,2,3])),2));
                    
                    ddddJ_Tdphidphidphidphi = ...
                        permute(temp,[1,3,2,4]) + permute(temp,[1,3,4,2]) + permute(temp,[3,4,1,2]);
                    
                    ddddJ_TdphidTdTdT = bsxfun(@times,ddddJ_TdTdTdTdT,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_TdTdTdTdSigma,permute(dSigma_timedphi,[3,1,2]));
                    ddddJ_TdphidTdTdSigma = bsxfun(@times,ddddJ_TdTdTdTdSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_TdTdTdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    ddddJ_TdphidTdSigmadSigma = bsxfun(@times,ddddJ_TdTdTdSigmadSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_TdTdSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    ddddJ_TdphidSigmadSigmadSigma = bsxfun(@times,ddddJ_TdTdSigmadSigmadSigma,permute(dYdphi,[3,1,2])) ...
                        + bsxfun(@times,ddddJ_TdSigmadSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    ddddJ_TdphidphidTdT = bsxfun(@times,ddddJ_TdphidTdTdT,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddddJ_TdphidTdTdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                    ddddJ_TdphidphidTdSigma = bsxfun(@times,ddddJ_TdphidTdTdSigma,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddddJ_TdphidTdSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                    ddddJ_TdphidphidSigmadSigma = bsxfun(@times,ddddJ_TdphidTdSigmadSigma,permute(dYdphi,[3,1,4,2])) ...
                        + bsxfun(@times,ddddJ_TdphidSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                    ddddJ_TdphidphidphidT = bsxfun(@times,ddddJ_TdphidphidTdT,permute(dYdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,ddddJ_TdphidphidTdSigma,permute(dSigma_timedphi,[3,1,4,5,2]));
                    ddddJ_TdphidphidphidSigma = bsxfun(@times,ddddJ_TdphidphidTdSigma,permute(dYdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,ddddJ_TdphidphidSigmadSigma,permute(dSigma_timedphi,[3,1,4,5,2]));
                    ddddJ_Tdphidphidphidphi = ddddJ_Tdphidphidphidphi + squeeze(sum(bsxfun(@times,ddddJ_TdphidphidphidT,permute(dYdphi,[3,1,4,5,6,2])) ...
                        + bsxfun(@times,ddddJ_TdphidphidphidSigma,permute(dSigma_timedphi,[3,1,4,5,6,2])),2));
                    
                    ddddJ_Tdbdphidphidphi = permute(sum(bsxfun(@times,ddddJ_Tdphidphidphidphi,permute(dphidb,[1,3,4,5,2])),1),[5,2,3,4,1]);
                    ddddJ_Tdbdbdphidphi = permute(sum(bsxfun(@times,ddddJ_Tdbdphidphidphi,permute(dphidb,[3,1,4,5,2])),2),[1,5,3,4,2]) ...
                        + permute(sum(bsxfun(@times,dddJ_Tdphidphidphi,permute(ddphidbdb,[1,4,5,2,3])),1),[4,5,2,3,1]);
                    ddddJ_Tdbdbdbdphi = permute(sum(bsxfun(@times,ddddJ_Tdbdbdphidphi,permute(dphidb,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                        + permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(ddphidbdb,[4,1,5,2,3])),2),[1,4,5,3,2]);
                    ddddJ_Tdbdbdbdb = permute(sum(bsxfun(@times,ddddJ_Tdbdbdbdphi,permute(dphidb,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                        + permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(ddphidbdb,[4,5,1,2,3])),3),[1,2,4,5,3]);
                    
                    ddddJdbdbdbdb = ddddJ_Ddbdbdbdb + ddddJ_Tdbdbdbdb + squeeze(ddddJ_bdbdbdbdb);
                    
                    varargout{14} = ddddJdbdbdbdb;
                    
                    %% ddddJdbdbdbdbeta
                    
                    ddddJ_Ddbdbdbdbeta = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbdphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                        + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                    
                    ddddJ_Tdbdbdbdbeta = permute(sum(bsxfun(@times,ddddJ_Tdbdbdbdphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                        + permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(ddphidbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                    
                    ddddJdbdbdbdbeta = ddddJ_Ddbdbdbdbeta + ddddJ_Ddbdbdbdbeta;
                    
                    varargout{15} = ddddJdbdbdbdbeta;
                    
                    %% ddddJdbdbdbetadbeta
                    
                    ddddJ_Ddbdbdbetadphi = permute(sum(bsxfun(@times,ddddJ_Ddbdbdphidphi,permute(dphidbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                        + permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(ddphidbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
                    ddddJ_Ddbdbdbetadbeta = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbetadphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                        + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                    
                    ddddJ_Tdbdbdbetadphi = permute(sum(bsxfun(@times,ddddJ_Tdbdbdphidphi,permute(dphidbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                        + permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(ddphidbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
                    ddddJ_Tdbdbdbetadbeta = permute(sum(bsxfun(@times,ddddJ_tdbdbdbetadphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                        + permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(ddphidbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                    
                    ddddJdbdbdbetadbeta = ddddJ_Ddbdbdbetadbeta + ddddJ_Tdbdbdbetadbeta;
                    
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

function varargout = simulateForSP(model,tout,phi,kappa)

% Simulate model
if(nargout<2)
    options_simu.sens = 0;
    sol = model(tout,phi,kappa,options_simu.sens);
    varargout{1} = sol.y;
else
    options_simu.sens = 1;
    sol = model(tout,phi,kappa,options_simu.sens);
    varargout{1} = sol.y;
    varargout{2} = sol.sy;
end
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
                    varargout{15} = transpose(60*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^6)) - 6./(Sigma(ind).^4));
                end
            end
        end
    end
end
end

function varargout = lognormal_noise(Y,Ym,Sigma,ind)
if nargout >=1
    % J_D
    varargout{1} = transpose(sum(0.5*((log(Y(ind)) - Ym(ind))./Sigma).^2 + 0.5*log(sqrt(2*pi)*(Sigma(ind).^2).*(Y(ind).^2))));
    if nargout >= 3
        % dJ_DdY
        varargout{2} = transpose((log(Y(ind)) - Ym(ind))./(Sigma(ind).^2).*(1./Y) + 1./Y(ind));
        % dJ_DdSigma
        varargout{3} = transpose(- (((log(Y(ind)) - Ym(ind)).^2)./(Sigma(ind).^3)) + 1./Sigma(ind));
        if nargout >= 4
            %ddJ_DdYdY
            varargout{4} = transpose((1-(log((ind)) - Ym(ind)))./(Sigma(ind).^2.*Y.^2) - 1./(Y(ind).^2));
            %ddJ_DdYdSigma
            varargout{5} = transpose(-2*(Y(ind) - Ym(ind))./(Sigma(ind).^3.*Y));
            %ddJ_DdSigmadSigma
            varargout{6} = transpose(3*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^4)) - 1./((ind)^2));
            if nargout >= 7
                %dddJ_DdYdYdY
                varargout{7} = transpose(- 3./(Sigma(ind).^2.*Y.^3)-2./(Y.^3));
                %dddJ_DdYdYdSigma
                varargout{8} = transpose((1-(log(Y(ind)) - Ym(ind)))./(Sigma(ind).^3.*Y.^2));
                %dddJ_DdYdSigmadSigma
                varargout{9} = transpose(+ 6*(Y(ind) - Ym(ind))./(Sigma(ind).^4.*Y));
                %dddJ_DdSigmadSigmadSigma
                varargout{10} = transpose(- 12*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^5)) + 2./(Sigma(ind)^3));
                if nargout >= 11
                    %ddddJ_DdYdYdYdY
                    varargout{11} = transpose(+9./(Sigma(ind).^2.*Y.^4)+6./(Y.^4));
                    %ddddJ_DdYdYdYdSigma
                    varargout{12} = transpose(6./(Sigma(ind).^3.*Y.^3)-2./(Y.^3));
                    %ddddJ_DdYdYdSigmadSigma
                    varargout{13} = transpose(-3*(1-(log(Y(ind)) - Ym(ind)))./(Sigma(ind).^4.*Y.^2));
                    %ddddJ_DdYdSigmadSigmadSigma
                    varargout{14} = transpose(- 24*(Y(ind) - Ym(ind))./(Sigma(ind).^5.*Y));
                    %ddddJ_DdSigmadSigmadSigmadSigma
                    varargout{15} = transpose(- 60*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^6)) + 6./(Sigma(ind)^4));
                end
            end
        end
    end
end
end

function varargout = normal_time(T,Tm,R,Sigma,ind)
if nargout >=1
    % J_D
    varargout{1} = sum(0.5*((T(ind) - Tm(ind))./Sigma(ind)).^2 + sum(0.5*((R(ind))./Sigma(ind)).^2) + log(sqrt(2*pi)*Sigma(ind).^2));
    if nargout >= 2
        % dJ_DdT
        varargout{2} = transpose((T(ind) - Tm(ind))./(Sigma(ind).^2) + R(ind)./(Sigma(ind).^2));
        % dJ_DdSigma
        varargout{3} = transpose(- (((T(ind) - Tm(ind)).^2)./(Sigma(ind).^3)) - (((R(ind)).^2)./(Sigma(ind).^3)) + 2./Sigma(ind));
        if nargout >= 4
            %ddJ_DdTdT
            varargout{4} = transpose(1./(Sigma(ind).^2));
            %ddJ_DdTdSigma
            varargout{5} = transpose(-2*(T(ind) - Tm(ind))./(Sigma(ind).^3) -2*(R(ind))./(Sigma(ind).^3));
            %ddJ_DdSigmadSigma
            varargout{6} = transpose(3*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^4)) + 3*(((R(ind)).^2)./(Sigma(ind).^4)) - 2./(Sigma(ind).^2));
            if nargout >= 7
                %dddJ_DdTdTdT
                varargout{7} = transpose(zeros(size(T(ind))));
                %dddJ_DdTdTdSigma
                varargout{8} = transpose(- 2./(Sigma(ind).^3));
                %dddJ_DdTdSigmadSigma
                varargout{9} = transpose(6*(T(ind) - Tm(ind))./(Sigma(ind).^4) + 6*(R(ind))./(Sigma(ind).^4));
                %dddJ_DdSigmadSigmadSigma
                varargout{10} = transpose(- 12*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^5)) - 12*(((R(ind)).^2)./(Sigma(ind).^5)) + 4./(Sigma(ind).^3));
                if nargout >= 11
                    %ddddJ_DdTdTdTdT
                    varargout{11} = transpose(zeros(size(T(ind))));
                    %ddddJ_DdTdTdTdSigma
                    varargout{12} = transpose(zeros(size(T(ind))));
                    %ddddJ_DdTdTdSigmadSigma
                    varargout{13} = transpose(6./(Sigma(ind).^4));
                    %ddddJ_DdTdSigmadSigmadSigma
                    varargout{14} = transpose(- 24*((T(ind) - Tm(ind))./(Sigma(ind).^5)) - 24*((R(ind))./(Sigma(ind).^5)));
                    %ddddJ_DdSigmadSigmadSigmadSigma
                    varargout{15} = transpose(60*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^6)) + 60*(((R(ind)).^2)./(Sigma(ind).^6)) - 12./(Sigma(ind).^4));
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




