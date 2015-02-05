% function [logL,dlogLdtheta,ddlogLdtheta2] = logL_CE_w_grad_2(theta,Data,Model,options)
function varargout = logL_CE_w_grad_3(varargin)
    
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
    theta = varargin{1};
    Data = varargin{2};
    Model = varargin{3};
    
    % Options
    options.tau_update = 0;
    options.plot = 1;
    if nargin >= 4
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
    dlogLdtheta = zeros(length(theta),1);
    ddlogLdtheta2 = zeros(length(theta));
    
    % Data types
    % - Single Cell Time-Lapse (SCTL)
    % - Population Average (PA)
    % - Single cell SnapsHot (SCSH)
    data_type = {'SCTL','PA','SCSH'};
    
    % Loop: Experimental conditions
    for s = 1:length(Data)
        
        %% Assignment of global variables
        A = Model.exp{s}.A;
        B = Model.exp{s}.B;
        ind_beta = Model.exp{s}.ind_beta;
        ind_D = Model.exp{s}.ind_D;
        n_beta = length(ind_beta);
        n_D = length(ind_D);
        n_b = size(B,2);
        type_D = Model.type_D;
        
        %% Construct fixed effects and covariance matrix
        beta = theta(ind_beta);
        dbeta = eye(length(theta)); dbeta = dbeta(ind_beta,:);
        [D,invD,dD,dinvD,~,HinvD] = xi2D(theta(ind_D),type_D);
        dD_full = zeros(size(dD,1),size(dD,1),length(theta));
        dD_full(:,:,ind_D) = dD;
        
        %% Construction of time vector
        t_s = [];
        for dtype = 1:length(data_type)
            if isfield(Data{s},data_type{dtype})
                t_s = union(eval(['Data{s}.' data_type{dtype} '.time']),t_s);
            end
        end
        
        %% Single cell time-lapse data - Individuals
        if isfield(Data{s},'SCTL')
            % Evaluation of time index set
            [~,ind_t] = ismember(Data{s}.SCTL.time,t_s);
            
            % Initialization
            Sim_SCTL = nan(size(Data{s}.SCTL.Y));
            
            % Loop: Indiviudal cells
            for i = 1:size(Data{s}.SCTL.Y,3)
                % Load single-cell data
                Ym_si = Data{s}.SCTL.Y(:,:,i);
                Sigma_si = Data{s}.SCTL.Sigma_Y(:,:,i);
                ind = find(~isnan(Ym_si));
                
                % Optimization of random effects
                options_fmincon = optimset('algorithm','trust-region-reflective',...
                    'display','off','GradObj','on',...
                    'MaxIter',100,... % 1000
                    'TolFun',1e-8,...
                    'TolFun',1e-8,...
                    'Hessian','user-supplied');
                if logL_old == -inf
                    bhat_si0 = zeros(size(B,2),1);
                else
                    bhat_si0 = P_old{s}.SCTL.bhat(:,i);
                end
                bhat_si0 = zeros(size(B,2),1);
                dbhat_si = zeros(n_b,n_b+n_D);
                [bhat_si,~,~,~,~,grad,hessian] = fmincon(...
                    @(b) objective_SCTL_s1(Model.exp{s}.model,beta,b,Data{s}.condition,invD,A,B,t_s,Ym_si,Sigma_si,ind),...
                    bhat_si0,[],[],[],[],-5*ones(n_b,1),5*ones(n_b,1),[],options_fmincon);
                P{s}.SCTL.bhat(:,i) = bhat_si;
                G = hessian;
                invG = pinv(hessian);
                
                % Gradient of optimal point
                [J,dJdtheta,F,dFdb,dFdbeta] = objective_SCTL_s1_full(Model.exp{s}.model,beta,bhat_si,Data{s}.condition,invD,dinvD,HinvD,A,B,t_s,Ym_si,Sigma_si,ind);
                
                %
                % derivative test
                %
                %[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(bhat_si,@(bhat_si) objective_SCTL_s1_full(Model.exp{s}.model,beta,bhat_si,Data{s}.condition,invD,dinvD,HinvD,A,B,t_s,Ym_si,Sigma_si,ind),1e-4,3,4)
                %[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(beta,@(beta) objective_SCTL_s1_full(Model.exp{s}.model,beta,bhat_si,Data{s}.condition,invD,dinvD,HinvD,A,B,t_s,Ym_si,Sigma_si,ind),1e-4,3,5)
                
                
                for l = 1:(n_beta+n_D)
                    W = [F([1:n_b,n_b+l],[1:n_b,n_b+l]),[zeros(n_b,1);1];[zeros(n_b,1);1]',0];
                    g = [zeros(n_b+1,1);1];
                    d = W\g;
                    dbhat_si(:,l) = d(1:n_b);
                end
                
                
                dGdbeta = zeros(size(B,1),size(B,1),n_beta);
                dGdD = zeros(size(B,1),size(B,1),n_D);
                
                for l = 1:n_beta
                    dGdbeta(1:n_b,1:n_b,l) = sum(bsxfun(@times,dFdb(1:n_b,1:n_b,:),permute(dbhat_si(:,l),[3,2,1])),3) + dFdbeta(1:n_b,1:n_b,l);
                end
                
                for l = 1:n_D
                    dGdD(1:n_b,1:n_b,l) = sum(bsxfun(@times,dFdb(1:n_b,1:n_b,:),permute(dbhat_si(:,n_beta+l),[3,2,1])),3) -invD*dD(:,:,l)*invD;
                end
                
                
                % Construct single-cell parameter
                phi_si = A*beta + B*bhat_si;
                
                % Simulate model
                [~,~,~,y,~,sy] = Model.exp{s}.model(t_s,phi_si,Data{s}.condition);
                
                % Evaluation of likelihood and likelihood gradient
                Y_si = y(ind_t,:);
                logL = logL ...
                    - 0.5*sum(((Ym_si(ind) - Y_si(ind))./Sigma_si(ind)).^2) ... % - g
                    - 0.5*log(det(D)) - 0.5*bhat_si'*invD*bhat_si ... % - log(det(D)) - p(b)
                    - 0.5*log(det(G)); % - log(det(G))
                
                if nargout >= 2
                    % beta
                    for k = 1:length(ind_beta)
                        sY_si_k = sum(bsxfun(@times,sy(ind_t,:,:),permute(A(:,k) + B*dbhat_si(:,k),[3,2,1])),3); %dydbeta_k = dydphi*dphidbeta_
                        
                        dlogLdtheta(ind_beta(k)) = dlogLdtheta(ind_beta(k)) ...
                            + sum(((Ym_si(ind) - Y_si(ind))./Sigma_si(ind).^2).*sY_si_k(ind)) ...  % dgdy*dydbeta_k
                            - dbhat_si(:,k)'*invD*bhat_si ... % - dp(b)dbeta
                            - 0.5*trace(invG*dGdbeta(:,:,k)); % dlog(det(G))dbeta
                    end
                    % D
                    for k = 1:length(ind_D)
                        sY_si_k = sum(bsxfun(@times,sy(ind_t,:,:),permute(B*dbhat_si(:,n_beta+k),[3,2,1])),3);
                        
                        dlogLdtheta(ind_D(k)) = dlogLdtheta(ind_D(k)) ...
                            + sum(((Ym_si(ind) - Y_si(ind))./Sigma_si(ind).^2).*sY_si_k(ind)) ... %dgdy*dydD_k
                            - 0.5*trace(invD*dD(:,:,k)) ... %dlog(det(D))dD
                            - 0.5*bhat_si'*dinvD(:,:,k)*bhat_si ... %dp(n)dD
                            - dbhat_si(:,n_beta+k)'*invD*bhat_si ... %dp(n)dD
                            - 0.5*trace(invG*dGdD(:,:,k)); % dlog(det(G))dD
                    end
                end
                
                % Assignment
                Sim_SCTL(:,:,i) = y(ind_t,:);
            end
            
            % Visulization
            if options.plot
                Model.exp{s}.plot(Data{s},Sim_SCTL,Model.exp{s}.fh);
            end
        end
        
        %% Single cell time-lapse data - Statistics
        if isfield(Data{s},'SCTLstat')
            % Simulation using sigma points
            [~,~,~,mz_SP,Cz_SP,~,~,~,~,~,dmz_SP,dCz_SP,~,~] = ...
                getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition),...
                Model.exp{s}.A,Model.exp{s}.B,beta,D,dbeta,dD_full);
            
            % Evaluation of likelihood, likelihood gradient and hessian
            % Mean
            logL_mz = - 0.5*sum(sum(((Data{s}.SCTLstat.mz - mz_SP)./Data{s}.SCTLstat.Sigma_mz).^2,1),2);
            dlogL_mzdtheta = sum(bsxfun(@times,(Data{s}.SCTLstat.mz - mz_SP)./Data{s}.SCTLstat.Sigma_mz.^2,dmz_SP),1)';
            wdmz_SP = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_mz,dmz_SP);
            ddlogL_mzdtheta2 = -wdmz_SP'*wdmz_SP;
            
            % Covariance
            logL_Cz = - 0.5*sum(sum(sum(((Data{s}.SCTLstat.Cz - Cz_SP)./Data{s}.SCTLstat.Sigma_Cz).^2,1),2),3);
            dlogL_Czdtheta = squeeze(sum(sum(bsxfun(@times,(Data{s}.SCTLstat.Cz - Cz_SP)./Data{s}.SCTLstat.Sigma_Cz.^2,dCz_SP),1),2));
            wdCz_SP = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_Cz,dCz_SP);
            wdCz_SP = reshape(wdCz_SP,[prod(size(Cz_SP)),size(dCz_SP,3)]);
            ddlogL_Czdtheta2 = -wdCz_SP'*wdCz_SP;
            
            % Summation
            logL = logL + logL_mz + logL_Cz;
            dlogLdtheta = dlogLdtheta + dlogL_mzdtheta + dlogL_Czdtheta;
            ddlogLdtheta2 = ddlogLdtheta2 + ddlogL_mzdtheta2 + ddlogL_Czdtheta2;
            
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
            [m_SP,C_SP,~,~,~,~,~,dm_SP,dC_SP,~,~,~,~,~] = ...
                getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCSH.time,phi,Data{s}.condition),...
                Model.exp{s}.A,Model.exp{s}.B,beta,D,dbeta,dD_full);
            
            % Evaluation of likelihood, likelihood gradient and hessian
            % Mean
            logL_m = - 0.5*sum(sum(((Data{s}.SCSH.m - m_SP)./Data{s}.SCSH.Sigma_m).^2,1),2);
            dlogL_mdtheta = squeeze(sum(sum(bsxfun(@times,(Data{s}.SCSH.m - m_SP)./Data{s}.SCSH.Sigma_m.^2,dm_SP),1),2));
            wdm_SP = bsxfun(@times,1./Data{s}.SCSH.Sigma_m,dm_SP);
            wdm_SP = reshape(wdm_SP,[prod(size(m_SP)),size(dm_SP,3)]);
            ddlogL_mdtheta2 = -wdm_SP'*wdm_SP;
            
            % Covariance
            logL_C = - 0.5*sum(sum(sum(((Data{s}.SCSH.C - C_SP)./Data{s}.SCSH.Sigma_C).^2,1),2),3);
            dlogL_Cdtheta = squeeze(sum(sum(sum(bsxfun(@times,(Data{s}.SCSH.C - C_SP)./Data{s}.SCSH.Sigma_C.^2,dC_SP),1),2),3));
            wdC_SP = bsxfun(@times,1./Data{s}.SCSH.Sigma_C,dC_SP);
            wdC_SP = reshape(wdC_SP,[prod(size(C_SP)),size(dC_SP,4)]);
            ddlogL_Cdtheta2 = -wdC_SP'*wdC_SP;
            
            % Summation
            logL = logL + logL_m + logL_C;
            dlogLdtheta = dlogLdtheta + dlogL_mdtheta + dlogL_Cdtheta;
            ddlogLdtheta2 = ddlogLdtheta2 + ddlogL_mdtheta2 + ddlogL_Cdtheta2;
            
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
            [m_SP,~,~,~,~,~,~,dm_SP,~,~,~,~,~,~] = ...
                getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.PA.time,phi,Data{s}.condition),...
                Model.exp{s}.A,Model.exp{s}.B,beta,D,dbeta,dD_full);
            
            % Evaluation of likelihood, likelihood gradient and hessian
            logL_m = - 0.5*sum(sum(((Data{s}.PA.m - m_SP)./Data{s}.PA.Sigma_m).^2,1),2);
            dlogL_mdtheta = squeeze(sum(sum(bsxfun(@times,(Data{s}.PA.m - m_SP)./Data{s}.PA.Sigma_m.^2,dm_SP),1),2));
            wdm_SP = bsxfun(@times,1./Data{s}.PA.Sigma_m,dm_SP);
            wdm_SP = reshape(wdm_SP,[prod(size(m_SP)),size(dm_SP,3)]);
            ddlogL_mdtheta2 = -wdm_SP'*wdm_SP;
            
            % Summation
            logL = logL + logL_m;
            dlogLdtheta = dlogLdtheta + dlogL_mdtheta;
            ddlogLdtheta2 = ddlogLdtheta2 + ddlogL_mdtheta2;
            
            % Visulization
            if options.plot
                Sim_PA.m = m_SP;
                Model.exp{s}.plot(Data{s},Sim_PA,Model.exp{s}.fh);
            end
        end
        
    end
    
    if extract_flag
        if isfield(Data{s},'SCTL')
            varargout{1} = P;
        elseif isfield(Data{s},'SCTLstat')
            varargout{1} = B_SP;
        end
        return
    end
    
    %% Output
    if nargout <= 1
        % One output
        varargout{1} =  logL;
    elseif nargout <= 2
        % Two outputs
        varargout{1} =  logL;
        varargout{2} =  dlogLdtheta;
    else
        % Two outputs
        varargout{1} =  logL;
        varargout{2} =  dlogLdtheta;
        varargout{3} =  ddlogLdtheta2;
    end
    
end

%%
function [J,gradJ,F] = objective_SCTL_s1(model,beta,b,kappa,invD,A,B,t,Ym,Sigma,ind)
    
    % Single-cell parameters
    phi = A*beta + B*b;
    
    % Simulate model
    [~,~,~,Y,~,sY] = model(t,phi,kappa);
    
    % Residual and residual gradient
    res = (Y(ind)-Ym(ind))./Sigma(ind);
    dres = zeros(length(ind),length(b));
    for k = 1:length(b)
        sY_k = sum(bsxfun(@times,sY,permute(B(:,k),[3,2,1])),3);
        dres(:,k) = sY_k(ind)./Sigma(ind);
    end
    
    % Evaluation of likelihood
    J = 0.5*res'*res + 0.5*b'*invD*b;
    gradJ = dres'*res + invD*b;
    F = dres'*dres + invD;
    
end

function [J,gradJ,F,dFdb,dFdbeta] = objective_SCTL_s1_full(model,beta,b,kappa,invD,dinvD,HinvD,A,B,t,Ym,Sigma,ind)
    
    % Dimensionality
    n_beta = length(beta);
    n_b = length(b);
    n_D = size(dinvD,3);
    
    % Single-cell parameters
    phi = A*beta + B*b;
    % beta_i = beta+[0;1e-3];
    % phi_i = A*beta_i + B*b;
    
    % Simulate model
    [~,~,~,Y,~,sY,~,s2Y] = model(t,phi,kappa);
    % [~,~,~,Y_i,~,sY_i,~,s2Y_i] = model(t,phi_i,kappa);
    % test derivatives
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi) model(t,phi,kappa),1e-4,4,6)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi) model(t,phi,kappa),1e-4,6,8)
    
    % Residual and residual gradient
    res = (Y(ind)-Ym(ind))./Sigma(ind);
    % res_i = (Y_i(ind)-Ym(ind))./Sigma(ind);
    % dres_i = zeros(length(ind),n_b+n_beta+n_D);
    % for k = 1:n_b
    %     sY_k_i = sum(bsxfun(@times,sY_i,permute(B(:,k),[3,2,1])),3);
    %     dres_i(:,k) = sY_k_i(ind)./Sigma(ind);
    % end
    % for k = 1:n_beta
    %     sY_k_i = sum(bsxfun(@times,sY_i,permute(A(:,k),[3,2,1])),3); % dydbeta_k = <dydphi,dphidbeta_k>
    %     dres_i(:,n_b+k) = sY_k_i(ind)./Sigma(ind); %dresdbeta
    % end
    dres = zeros(length(ind),n_b+n_beta+n_D);
    d2res = zeros(length(ind),n_b+n_beta+n_D,n_b+n_beta+n_D);
    
    
    for k = 1:n_b
        sY_k = sum(bsxfun(@times,sY,permute(B(:,k),[3,2,1])),3); % dydb_k = sum_j d{y}d{phi_j} * d{phi_j}d{b_k}>
        dres(:,k) = sY_k(ind)./Sigma(ind); %dresdb
        for l = 1:k
            s2Y_kl = sum(sum(bsxfun(@times,s2Y,permute(B(:,k)*B(:,l)',[4,3,1,2])),4),3); %d2{y}d{b_k}d{b_l} = sum_j d2{y}d{phi_k}d{phi_j} * d{phi_j}d{b_l}
            s2Y_lk = sum(sum(bsxfun(@times,s2Y,permute(B(:,l)*B(:,k)',[4,3,1,2])),4),3);
            d2res(:,k,l) = s2Y_kl(ind)./Sigma(ind);
            d2res(:,l,k) = s2Y_lk(ind)./Sigma(ind);
        end
        for l = 1:n_beta
            s2Y_kl = sum(sum(bsxfun(@times,s2Y,permute(B(:,k)*A(:,l)',[4,3,1,2])),4),3); %d2{y}d{b_k}d{beta_l} = sum_j d2{y}d{phi_k}d{phi_j} * d{phi_j}d{b_l}
            s2Y_lk = sum(sum(bsxfun(@times,s2Y,permute(A(:,l)*B(:,k)',[4,3,1,2])),4),3);
            d2res(:,k,n_b+l) = s2Y_kl(ind)./Sigma(ind);
            d2res(:,n_b+l,k) = s2Y_lk(ind)./Sigma(ind);
        end
    end
    for k = 1:n_beta
        sY_k = sum(bsxfun(@times,sY,permute(A(:,k),[3,2,1])),3); % dydbeta_k = <dydphi,dphidbeta_k>
        dres(:,n_b+k) = sY_k(ind)./Sigma(ind); %dresdbeta
        for l = 1:n_beta
            s2Y_kl = sum(sum(bsxfun(@times,s2Y,permute(A(:,k)*A(:,l)',[4,3,1,2])),4),3); %d2{y}d{beta_k}d{beta_l} = sum_j sum_m d2{y}d{phi_j}d{phi_m} * d{phi_j}d{beta_k} * d{phi_m}d{beta_l}
            s2Y_lk = sum(sum(bsxfun(@times,s2Y,permute(A(:,l)*A(:,k)',[4,3,1,2])),4),3);
            d2res(:,n_b+k,n_b+l) = s2Y_kl(ind)./Sigma(ind);
            d2res(:,n_b+l,n_b+k) = s2Y_lk(ind)./Sigma(ind);
        end
    end
    
    % Covariance components
    grad_D_component = zeros(n_D,1);
    H_D_component = zeros(n_D,n_D);
    H_Db_component = zeros(n_b,n_D);
    for k1 = 1:n_D
        grad_D_component(k1) = 0.5*b'*dinvD(:,:,k1)*b;
        H_Db_component(:,k1) = dinvD(:,:,k1)*b;
        for k2 = 1:k1
            H_D_component(k1,k2) = b'*HinvD(:,:,k1,k2)*b;
            H_D_component(k2,k1) = b'*HinvD(:,:,k2,k1)*b;
        end
    end
    
    
    % Evaluation of likelihood
    J = 0.5*res'*res + 0.5*b'*invD*b;
    gradJ = dres'*res+[invD*b;...
        zeros(n_beta,1);...
        grad_D_component];
    F = dres'*dres + ...
        squeeze(sum(bsxfun(@times,res,d2res),1)) + ...
        [invD,zeros(n_b,n_beta),H_Db_component;...
        zeros(n_beta,n_b+n_beta+n_D);...
        H_Db_component',zeros(n_D,n_beta),H_D_component];
    
    
    dFdb = zeros(n_b + n_beta + n_D,n_b + n_beta + n_D,n_b);
    dFdbeta = zeros(n_b + n_beta + n_D,n_b + n_beta + n_D,n_beta);
    
    for l = 1:n_b + n_beta + n_D
        for k = 1:1:n_b + n_beta + n_D
            for i = 1:n_b
                dFdb(l,k,i) = sum(bsxfun(@times,d2res(:,l,i),dres(:,k)),1) + ...
                    sum(bsxfun(@times,d2res(:,k,l),dres(:,i)),1) + ...
                    sum(bsxfun(@times,d2res(:,i,k),dres(:,l)),1);
            end
            for i = 1:n_beta
                dFdbeta(l,k,i) = sum(bsxfun(@times,d2res(:,l,n_b+i),dres(:,k)),1) + ...
                    sum(bsxfun(@times,d2res(:,k,l),dres(:,n_b+i)),1) + ...
                    sum(bsxfun(@times,d2res(:,n_b+i,k),dres(:,l)),1);
            end
        end
    end
    
    % H_D_component
    for l = 1:n_D
        for k = 1:n_D
            for i = 1:n_b
                dFdb(n_b+n_beta+l,n_b+n_beta+k,i) = 2*sum(HinvD(:,:,k,l)*b,1);
            end
        end
    end
    % H_Db_component
    for l = 1:n_D
        for k = 1:n_b
            for i = 1:n_b
                dFdb(l,n_b+n_beta+k,i) = dinvD(k,i,l);
                dFdb(n_b+n_beta+l,k,i) = dinvD(i,k,l);
            end
        end
    end
    
end

%%
function [y,sy] = simulateForSP(model,tout,phi,kappa)
    
    % Simulate model
    [status,t,x,y,sx,sy] = model(tout,phi,kappa);
    
end