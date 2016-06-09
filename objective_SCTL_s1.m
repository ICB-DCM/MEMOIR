% objective_SCTL_s1 is an auxiliary function for logL_CE_w_grad_2 and
% serves as objective function for the estimation of random effect
% parameters
%
% USAGE:
% ======
% [J,J.db,J.dbdb,...] = objective_SCTL_s1(model,beta,b,data.condition,delta,options.type_D,t,data.SCTL.Y(:,:,i),data.SCTL.T(:,:,i),data.SCTL.ind_y(:,i),data.SCTL.ind_t(:,i),s)
%
% INPUTS:
% =======
%
% model ... model definition
% beta ... common effect parameter
% b ... random effect parameter
% data.condition ... experimental condition
% delta ... parametrisation of random effect covariance
% options.type_D ... covariance parametrisation type
% t ... time vector for simulation
% data.SCTL.Y(:,:,i) ... measurements
% data.SCTL.T(:,:,i) ... observed event-times
% data.SCTL.ind_y(:,i) ... indexing of measurements
% data.SCTL.ind_t(:,i) ... indexing of events
% s ... experimental index
%
% Outputs:
% ========
% objective function J and derivatives wrt to b, beta and delta, see file
% for details, pd indicates that only the partial derivative is considered
% J
% J.db
% J.dbdb = G
% J.dbdbeta
% J.dbddelta
% J.dbetadbeta
% J.ddeltaddelta
% J.dbetaddelta
% dGdb
% pdGpdbeta
% pdGpddelta
% J.dbdbetadbeta
% J.dbddeltaddelta
% J.dbdbetaddelta
% ddGdbdb
% pddGdbpdbeta
% pdpdGpdbetapdbeta
% pddGdbpddelta
% pdpdGpddeltapddelta
% pdpdGpdbetapddelta
%
% 2015/04/14 Fabian Froehlich

function [varargout] = objective_SCTL_s1(model,data,beta,b,delta,s,i,options,nderiv)
    
    [D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,options.type_D);
    
    % mixed effect parameter
    phi = model.phi(beta,b);
    
    nderiv = max(nderiv,nargout>1);
    
    % build standard deviation matrices
    switch(nderiv)
        case 0
        [Sigma_noise] = build_sigma_noise(phi,data.SCTL.Y(:,:,i),s,model,data.SCTL.ind_y(:,i));
        [Sigma_time] = build_sigma_time(phi,data.SCTL.T(:,:,i),s,model,data.SCTL.ind_t(:,i));
        case 1
        [Sigma_noise,dSigma_noisedphi] = build_sigma_noise(phi,data.SCTL.Y(:,:,i),s,model,data.SCTL.ind_y(:,i));
        [Sigma_time,dSigma_timedphi] = build_sigma_time(phi,data.SCTL.T(:,:,i),s,model,data.SCTL.ind_t(:,i));
        case 2
        [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi] = build_sigma_noise(phi,data.SCTL.Y(:,:,i),s,model,data.SCTL.ind_y(:,i));
        [Sigma_time,dSigma_timedphi,ddSigma_timedphidphi] = build_sigma_time(phi,data.SCTL.T(:,:,i),s,model,data.SCTL.ind_t(:,i));
        case 3
        [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi] = build_sigma_noise(phi,data.SCTL.Y(:,:,i),s,model,data.SCTL.ind_y(:,i));
        [Sigma_time,dSigma_timedphi,ddSigma_timedphidphi,dddSigma_timedphidphidphi] = build_sigma_time(phi,data.SCTL.T(:,:,i),s,model,data.SCTL.ind_t(:,i));
        case 4
        [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi,ddddSigma_noisedphidphidphidphi] = build_sigma_noise(phi,data.SCTL.Y(:,:,i),s,model,data.SCTL.ind_y(:,i));
        [Sigma_time,dSigma_timedphi,ddSigma_timedphidphi,dddSigma_timedphidphidphi,ddddSigma_timedphidphidphidphi] = build_sigma_time(phi,data.SCTL.T(:,:,i),s,model,data.SCTL.ind_t(:,i));
    end
    
    t = data.SCTL.time(data.SCTL.ind_y(:,i));
    if(t(end)<max(max(data.SCTL.T(:,:,i)))*1.2)
        t = [t;max(max(data.SCTL.T(:,:,i)))*1.2];
    end
    
    % simulate trajectory
    if(nargout == 1)
        [Y,T,R] = simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i));
    elseif(or(nargout == 2,nargout == 3))
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i)),1e-5,1,4)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i)),1e-5,2,5)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i)),1e-5,3,6)
        [Y,T,R,dYdphi,dTdphi,dRdphi] = simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i));
    else
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i)),1e-5,4,7)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i)),1e-5,5,8)
        % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi)simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i)),1e-5,6,9)
        [Y,T,R,dYdphi,dTdphi,dRdphi,ddYdphidphi,ddTdphidphi,ddRdphidphi] = simulate_trajectory(t,phi,model,data.condition,s,data.SCTL.ind_t(:,i),data.SCTL.ind_y(:,i));
    end
    
    % noise model
    if(~isfield(model,'noise_model'))
        model.noise_model = 'normal';
    end
    
    switch(model.noise_model)
        case 'normal'
            noisedist = @normal_noise;
        case 'lognormal'
            noisedist = @lognormal_noise;
        case 'tdist'
            noisedist = @tdist_noise;
    end
    switch(nderiv)
        case 0
            J_D = noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i));
        case 1 % first order derivatives
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Y,@(Y)noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i)),1e-3,1,2)
            % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_noise,@(Sigma_noise)noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i)),1e-5,1,3)
            if(nargout<3)
                [J_D,...
                    dJ_DdY,dJ_DdSigma] = noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i));
            else
                [J_D,...
                    dJ_DdY,dJ_DdSigma,...
                    ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i));
            end
        case 2 % second order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma] = noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i));
        case 3 % third order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma] = noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i));
        case 4% fourth order derivatives
            [J_D,...
                dJ_DdY,dJ_DdSigma,...
                ddJ_DdYdY,ddJ_DdYdSigma,ddJ_DdSigmadSigma,...
                dddJ_DdYdYdY,dddJ_DdYdYdSigma,dddJ_DdYdSigmadSigma,dddJ_DdSigmadSigmadSigma,...
                ddddJ_DdYdYdYdY,ddddJ_DdYdYdYdSigma,ddddJ_DdYdYdSigmadSigma,ddddJ_DdYdSigmadSigmadSigma,ddddJ_DdSigmadSigmadSigmadSigma] = noisedist(Y,data.SCTL.Y(:,:,i),Sigma_noise,data.SCTL.ind_y(:,i));
    end
    
    % event model
    if(~isfield(model,'time_model'))
        model.time_model = 'normal';
    end
    
    switch(model.time_model)
        case 'normal'
            switch(nderiv)
                case 0
                    J_T = normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i));
                case 1 % first order derivatives
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(T,@(T)normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i)),1e-5,1,2)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(R,@(R)normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i)),1e-5,1,3)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(Sigma_time,@(Sigma_time)normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i)),1e-5,1,4)
                    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(T,@(T)normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i)),1e-5,2,5)
                    if nargout<3
                        [J_T,...
                            dJ_TdT,dJ_TdR,dJ_TdSigma] = normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i));
                    else
                        [J_T,...
                            dJ_TdT,dJ_TdR,dJ_TdSigma,...
                            ddJ_TdTdT,ddJ_TdTdR,ddJ_TdRdR,ddJ_TdTdSigma,ddJ_TdRdSigma,ddJ_TdSigmadSigma] = normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i));
                    end
                case 2 % second order derivatives
                    [J_T,...
                        dJ_TdT,dJ_TdR,dJ_TdSigma,...
                        ddJ_TdTdT,ddJ_TdTdR,ddJ_TdRdR,ddJ_TdTdSigma,ddJ_TdRdSigma,ddJ_TdSigmadSigma] = normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i));
                case 3 % third order derivatives
                    [J_T,...
                        dJ_TdT,dJ_TdR,dJ_TdSigma,...
                        ddJ_TdTdT,ddJ_TdTdR,ddJ_TdRdR,ddJ_TdTdSigma,ddJ_TdRdSigma,ddJ_TdSigmadSigma,...
                        dddJ_TdTdTdT,dddJ_TdTdTdR,dddJ_TdTdRdR,dddJ_TdRdRdR,dddJ_TdTdTdSigma,dddJ_TdTdRdSigma,dddJ_TdRdRdSigma,dddJ_TdTdSigmadSigma,dddJ_TdRdSigmadSigma,dddJ_TdSigmadSigmadSigma] = normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i));
                case 4 % fourth order derivatives
                    [J_T,...
                        dJ_TdT,dJ_TdR,dJ_TdSigma,...
                        ddJ_TdTdT,ddJ_TdTdR,ddJ_TdRdR,ddJ_TdTdSigma,ddJ_TdRdSigma,ddJ_TdSigmadSigma,...
                        dddJ_TdTdTdT,dddJ_TdTdTdR,dddJ_TdTdRdR,dddJ_TdRdRdR,dddJ_TdTdTdSigma,dddJ_TdTdRdSigma,dddJ_TdRdRdSigma,dddJ_TdTdSigmadSigma,dddJ_TdRdSigmadSigma,dddJ_TdSigmadSigmadSigma,...
                        ddddJ_TdTdTdTdT,ddddJ_TdTdTdTdR,ddddJ_TdTdTdRdR,ddddJ_TdTdRdRdR,ddddJ_TdRdRdRdR,ddddJ_TdTdTdTdSigma,ddddJ_TdTdTdRdSigma,ddddJ_TdTdRdRdSigma,ddddJ_TdRdRdRdSigma,ddddJ_TdTdTdSigmadSigma,ddddJ_TdTdRdSigmadSigma,ddddJ_TdRdRdSigmadSigma,ddddJ_TdTdSigmadSigmadSigma,ddddJ_TdRdSigmadSigmadSigma,ddddJ_TdSigmadSigmadSigmadSigma] = normal_time(T,data.SCTL.T(:,:,i),R,Sigma_time,data.SCTL.ind_t(:,i));
            end
            
    end
    
    % random effect model
    switch(model.parameter_model)
        case 'normal'
            switch(nderiv)
                case 0
                    J_b = normal_param(b,delta,options.type_D);
                case 1 % first order derivatives
                    if(nargout<3)
                        [J_b,dJ_bdb,dJ_bddelta]= normal_param(b,delta,options.type_D);
                    else
                        [J_b,dJ_bdb,dJ_bddelta,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta] = normal_param(b,delta,options.type_D);
                    end
                case 2% second order derivatives
                    [J_b,dJ_bdb,dJ_bddelta,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta] = normal_param(b,delta,options.type_D);
                case 3 % third order derivatives
                    [J_b,dJ_bdb,dJ_bddelta,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta,dddJ_bdbdbdb,dddJ_bdbdbddelta,dddJ_bdbddeltaddelta] = normal_param(b,delta,options.type_D);
                case 4 % fourth order derivatives
                    [J_b,dJ_bdb,dJ_bddelta,ddJ_bdbdb,ddJ_bdbddelta,ddJ_bddeltaddelta,dddJ_bdbdbdb,dddJ_bdbdbddelta,dddJ_bdbddeltaddelta,ddddJ_bdbdbdbdb,ddddJ_bdbdbdbddelta,ddddJ_bdbdbddeltaddelta] = normal_param(b,delta,options.type_D);
            end
        case 'lognormal'
    end
    
    J.val = J_D + J_T + J_b ;
    
    if nargout >= 2
        %% J.db
        dphidb = model.dphidb(beta,b);
        dphidbeta = model.dphidbeta(beta,b);
        
        dJ_Ddphi = chainrule(dJ_DdY,dYdphi) + chainrule(dJ_DdSigma,dSigma_noisedphi);
        
        dJ_Tdphi = chainrule(dJ_TdT,dTdphi) + chainrule(dJ_TdR,dRdphi) + chainrule(dJ_TdSigma,dSigma_timedphi);
        
        dJ_Ddb = chainrule(dJ_Ddphi,dphidb);
        
        dJ_Tdb = chainrule(dJ_Tdphi,dphidb);
        
        J.db = dJ_Ddb + dJ_Tdb + dJ_bdb;
        
        %% J.dbeta
        
        dJ_Ddbeta = chainrule(dJ_Ddphi,dphidbeta);
        
        dJ_Tdbeta = chainrule(dJ_Tdphi,dphidbeta);
        
        J.dbeta = dJ_Ddbeta + dJ_Tdbeta;
        
        %% J.ddelta
        
        J.ddelta = dJ_bddelta;
        
        if nargout >= 3 || nderiv >= 1
            %% J.dbdb
            % we need to make two different computations here,
            % one for the integration, in order to ensure that it is possible
            % to compute derivatives by using second order sensitivities by
            % using the FIM approximation and one that is used in the
            % computation of implicit derivatives where the accuracy is more
            % important
            
            ddJ_DdphidY = bsxfun(@times,ddJ_DdYdY,permute(dYdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_DdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
            ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dYdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
            % approximate
            ddJ_Dappdphidphi = permute(nansum(bsxfun(@times,ddJ_DdphidY,permute(dYdphi,[3,1,4,2])) ...
                + bsxfun(@times,ddJ_DdphidSigma,permute(dSigma_noisedphi,[3,1,4,2])),2),[3,4,1,2]);
            
            
            ddphidbdb = model.ddphidbdb(beta,b);
            
            ddJ_Dappdbdphi = transpose(squeeze(nansum(bsxfun(@times,ddJ_Dappdphidphi,permute(dphidb,[3,1,2])),2)));
            ddJ_Dappdbdb = squeeze(sum(bsxfun(@times,ddJ_Dappdbdphi,permute(dphidb,[3,1,2,4])),2)) ...
                + chainrule(dJ_Ddphi,ddphidbdb);
            
            ddJ_TdphidT = bsxfun(@times,ddJ_TdTdT,permute(dTdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_TdTdR,permute(dRdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_TdTdSigma,permute(dSigma_timedphi,[3,1,2]));
            ddJ_TdphidR = bsxfun(@times,ddJ_TdTdR,permute(dTdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_TdRdR,permute(dRdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_TdRdSigma,permute(dSigma_timedphi,[3,1,2]));
            ddJ_TdphidSigma = bsxfun(@times,ddJ_TdTdSigma,permute(dTdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_TdRdSigma,permute(dRdphi,[3,1,2])) ...
                + bsxfun(@times,ddJ_TdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
            ddJ_Tappdphidphi = permute(sum(bsxfun(@times,ddJ_TdphidT,permute(dTdphi,[3,1,4,2])) ...
                + bsxfun(@times,ddJ_TdphidR,permute(dRdphi,[3,1,4,2])) ...
                + bsxfun(@times,ddJ_TdphidSigma,permute(dSigma_timedphi,[3,1,4,2])),2),[3,4,1,2]);
            
            ddphidbdb = model.ddphidbdb(beta,b);
            
            ddJ_Tappdbdphi = transpose(squeeze(nansum(bsxfun(@times,ddJ_Tappdphidphi,permute(dphidb,[3,1,2])),2)));
            ddJ_Tappdbdb = squeeze(nansum(bsxfun(@times,ddJ_Tappdbdphi,permute(dphidb,[3,1,2,4])),2)) + chainrule(dJ_Tdphi,ddphidbdb);
            
            % FIM
            FIM.val = squeeze(ddJ_Dappdbdb) + squeeze(ddJ_Tappdbdb) + squeeze(ddJ_bdbdb);
            
            if nderiv >= 2
                
                % exact
                ddJ_Ddphidphi = ddJ_Dappdphidphi ...
                    + permute(nansum(bsxfun(@times,dJ_DdY,permute(ddYdphidphi,[4,1,2,3])),2),[3,4,1,2]);
                
                ddJ_Ddbdphi = transpose(squeeze(nansum(bsxfun(@times,ddJ_Ddphidphi,permute(dphidb,[3,1,2])),2)));
                ddJ_Ddbdb = squeeze(nansum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidb,[3,1,2,4])),2)) ...
                    + chainrule(dJ_Ddphi,ddphidbdb);
                
                ddJ_Tdphidphi = ddJ_Tappdphidphi ...
                    + permute(nansum(bsxfun(@times,dJ_TdT,permute(ddTdphidphi,[4,1,2,3])),2),[3,4,1,2]);
                
                ddJ_Tdbdphi = transpose(squeeze(nansum(bsxfun(@times,ddJ_Tdphidphi,permute(dphidb,[3,1,2])),2)));
                ddJ_Tdbdb = squeeze(nansum(bsxfun(@times,ddJ_Tdbdphi,permute(dphidb,[3,1,2,4])),2)) + chainrule(dJ_Tdphi,ddphidbdb);
                
                J.dbdb = squeeze(ddJ_Ddbdb) + squeeze(ddJ_Tdbdb) + squeeze(ddJ_bdbdb);
                
                %% J.dbdbeta
                ddphidbdbeta = model.ddphidbdbeta(beta,b);
                
                % if size of b == 1 then we have to permute second term,
                if(numel(b) == 1)
                    ddJ_Ddbdbeta = permute(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidbeta,[3,1,2])),2),[1,3,4,2]) + transpose(chainrule(dJ_Ddphi,ddphidbdbeta));
                    ddJ_Tdbdbeta = permute(sum(bsxfun(@times,ddJ_Tdbdphi,permute(dphidbeta,[3,1,2])),2),[1,3,4,2]) + transpose(chainrule(dJ_Tdphi,ddphidbdbeta));
                else
                    ddJ_Ddbdbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbdphi,permute(dphidbeta,[3,1,2])),2)) + chainrule(dJ_Ddphi,ddphidbdbeta);
                    ddJ_Tdbdbeta = squeeze(sum(bsxfun(@times,ddJ_Tdbdphi,permute(dphidbeta,[3,1,2])),2)) + chainrule(dJ_Tdphi,ddphidbdbeta);
                end
                
                J.dbdbeta = ddJ_Ddbdbeta + ddJ_Tdbdbeta;
                
                %% J.dbdelta
                
                J.dbddelta = ddJ_bdbddelta;
                
                %% J.dbetadbeta
                
                ddphidbetadbeta = model.ddphidbetadbeta(beta,b);
                
                ddJ_Ddbetadphi = permute(sum(bsxfun(@times,ddJ_Ddphidphi,permute(dphidbeta,[3,1,2])),2),[3,2,1]);
                ddJ_Tdbetadphi = permute(sum(bsxfun(@times,ddJ_Tdphidphi,permute(dphidbeta,[3,1,2])),2),[3,2,1]);
                if(numel(b)==1)
                    ddJ_Ddbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule(dJ_Ddphi,ddphidbetadbeta);
                    ddJ_Tdbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Tdbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule(dJ_Tdphi,ddphidbetadbeta);
                else
                    ddJ_Ddbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Ddbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule(dJ_Ddphi,ddphidbetadbeta);
                    ddJ_Tdbetadbeta = squeeze(sum(bsxfun(@times,ddJ_Tdbetadphi,permute(dphidbeta,[3,1,2])),2)) + chainrule(dJ_Tdphi,ddphidbetadbeta);
                end
                
                J.dbetadbeta = ddJ_Ddbetadbeta + ddJ_Tdbetadbeta;
                
                %% J.ddeltaddelta
                
                J.ddeltaddelta = ddJ_bddeltaddelta;
                
                %% J.dbetaddelta
                
                J.dbetaddelta = zeros(length(beta),size(dDddelta,3));
                
                if nderiv >= 3
                    
                    %% FIM.db
                    
                    temp = squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(ddYdphidphi,[4,1,5,2,3])) ...
                        + bsxfun(@times,ddJ_DdphidSigma,permute(ddSigma_noisedphidphi,[4,1,5,2,3])),2));
                    
                    dddJ_Ddphidphidphi = chainrule(dJ_DdSigma,dddSigma_noisedphidphidphi) ...
                        + permute(temp,[2,1,3]) + permute(temp,[1,2,3]);
                    
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
                        + squeeze(chainrule(ddJ_Ddbdphi,ddphidbdb));
                    
                    temp = squeeze(sum(bsxfun(@times,ddJ_TdphidT,permute(ddTdphidphi,[4,1,5,2,3])) ...
                        + bsxfun(@times,ddJ_TdphidR,permute(ddRdphidphi,[4,1,5,2,3])) ...
                        + bsxfun(@times,ddJ_TdphidSigma,permute(ddSigma_timedphidphi,[4,1,5,2,3])),2));
                    
                    dddJ_Tdphidphidphi = chainrule(dJ_TdSigma,dddSigma_timedphidphidphi) ...
                        + permute(temp,[2,1,3]) + permute(temp,[1,2,3]);
                    
                    dddJ_TdphidTdT = bsxfun(@times,dddJ_TdTdTdT,permute(dTdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdTdTdR,permute(dRdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdTdTdSigma,permute(dSigma_timedphi,[3,1,2]));
                    dddJ_TdphidTdR = bsxfun(@times,dddJ_TdTdTdR,permute(dTdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdTdRdR,permute(dRdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdTdRdSigma,permute(dSigma_timedphi,[3,1,2]));
                    dddJ_TdphidRdR = bsxfun(@times,dddJ_TdTdRdR,permute(dTdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdRdRdR,permute(dRdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdRdRdSigma,permute(dSigma_timedphi,[3,1,2]));
                    dddJ_TdphidTdSigma = bsxfun(@times,dddJ_TdTdTdSigma,permute(dTdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdTdRdSigma,permute(dRdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdTdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    dddJ_TdphidRdSigma = bsxfun(@times,dddJ_TdTdRdSigma,permute(dTdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdRdRdSigma,permute(dRdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdRdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    dddJ_TdphidSigmadSigma = bsxfun(@times,dddJ_TdTdSigmadSigma,permute(dTdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdRdSigmadSigma,permute(dRdphi,[3,1,2])) ...
                        + bsxfun(@times,dddJ_TdSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                    dddJ_TdphidphidT = bsxfun(@times,dddJ_TdphidTdT,permute(dTdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_TdphidTdR,permute(dRdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_TdphidTdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                    dddJ_TdphidphidR = bsxfun(@times,dddJ_TdphidTdR,permute(dTdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_TdphidRdR,permute(dRdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_TdphidRdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                    dddJ_TdphidphidSigma = bsxfun(@times,dddJ_TdphidTdSigma,permute(dTdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_TdphidRdSigma,permute(dRdphi,[3,1,4,2])) ...
                        + bsxfun(@times,dddJ_TdphidSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                    
                    dddJ_Tdphidphidphi = dddJ_Tdphidphidphi + squeeze(sum(bsxfun(@times,dddJ_TdphidphidT,permute(dTdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,dddJ_TdphidphidR,permute(dRdphi,[3,1,4,5,2])) ...
                        + bsxfun(@times,dddJ_TdphidphidSigma,permute(dSigma_timedphi,[3,1,4,5,2])),2));
                    
                    dddJ_Tdbdphidphi = permute(sum(bsxfun(@times,dddJ_Tdphidphidphi,permute(dphidb,[1,3,4,2])),1),[4,2,3,1]);
                    dddJ_Tdbdbdphi = permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(dphidb,[3,1,4,2])),2),[1,4,3,2]) ...
                        + permute(sum(bsxfun(@times,ddJ_Tdphidphi,permute(ddphidbdb,[1,4,2,3])),1),[3,4,2,1]);
                    dddJ_Tdbdbdb = permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(dphidb,[3,4,1,2])),3),[1,2,4,3]) ...
                        + squeeze(chainrule(ddJ_Tdbdphi,ddphidbdb));
                    
                    FIM.db = dddJ_Ddbdbdb + dddJ_Tdbdbdb + squeeze(dddJ_bdbdbdb);
                    
                    %% FIM.dbeta
                    
                    dddJ_Ddbdbdbeta = permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                    
                    dddJ_Tdbdbdbeta = permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                    
                    FIM.dbeta = dddJ_Ddbdbdbeta + dddJ_Tdbdbdbeta;
                    
                    %% FIM.ddelta
                    FIM.ddelta = dddJ_bdbdbddelta;
                    
                    %% J.dbdbetadbeta
                    
                    dddJ_Ddbdbetadphi = permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(dphidbeta,[3,1,4,2])),2),[1,4,3,2]);
                    dddJ_Ddbdbetadbeta = permute(sum(bsxfun(@times,dddJ_Ddbdbetadphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                    
                    dddJ_Tdbdbetadphi = permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(dphidbeta,[3,1,4,2])),2),[1,4,3,2]);
                    dddJ_Tdbdbetadbeta = permute(sum(bsxfun(@times,dddJ_Tdbdbetadphi,permute(dphidbeta,[3,4,1,2])),3),[1,2,4,3]);
                    
                    J.dbdbetadbeta = dddJ_Ddbdbetadbeta + dddJ_Tdbdbetadbeta;
                    
                    %% J.dbddeltadelta
                    J.dbddeltaddelta = dddJ_bdbddeltaddelta;
                    
                    %% J.dbdbetaddelta
                    
                    J.dbdbetaddelta = zeros(length(b),length(beta),size(dDddelta,3));
                    
                    if nderiv >= 4
                        
                        %% FIM.dbdb
                        
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
                        
                        temp = squeeze(sum(bsxfun(@times,bsxfun(@times,ddJ_TdTdT,permute(ddTdphidphi,[4,1,2,3])),permute(ddTdphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_TdTdR,permute(ddTdphidphi,[4,1,2,3])),permute(ddRdphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_TdRdR,permute(ddRdphidphi,[4,1,2,3])),permute(ddRdphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_TdTdSigma,permute(ddTdphidphi,[4,1,2,3])),permute(ddSigma_timedphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_TdRdSigma,permute(ddRdphidphi,[4,1,2,3])),permute(ddSigma_timedphidphi,[4,1,5,6,2,3])) ...
                            + bsxfun(@times,bsxfun(@times,ddJ_TdSigmadSigma,permute(ddSigma_timedphidphi,[4,1,2,3])),permute(ddSigma_timedphidphi,[4,1,5,6,2,3])),2));
                        
                        ddddJ_Tdphidphidphidphi = ...
                            permute(temp,[1,3,2,4]) + permute(temp,[1,3,4,2]) + permute(temp,[3,4,1,2]);
                        
                        ddddJ_TdphidTdTdT = bsxfun(@times,ddddJ_TdTdTdTdT,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdTdTdR,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdTdTdSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidTdTdR = bsxfun(@times,ddddJ_TdTdTdTdR,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdTdRdR,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdTdRdSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidTdRdR = bsxfun(@times,ddddJ_TdTdTdRdR,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdRdRdR,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdRdRdSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidRdRdR = bsxfun(@times,ddddJ_TdTdRdRdR,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdRdRdRdR,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdRdRdRdSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidTdTdSigma = bsxfun(@times,ddddJ_TdTdTdTdSigma,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdTdRdSigma,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdTdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidTdRdSigma = bsxfun(@times,ddddJ_TdTdTdRdSigma,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdRdRdSigma,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdRdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidRdRdSigma = bsxfun(@times,ddddJ_TdTdRdRdSigma,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdRdRdRdSigma,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdRdRdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidTdSigmadSigma = bsxfun(@times,ddddJ_TdTdTdSigmadSigma,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdRdSigmadSigma,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidRdSigmadSigma = bsxfun(@times,ddddJ_TdTdRdSigmadSigma,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdRdRdSigmadSigma,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdTdSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidSigmadSigmadSigma = bsxfun(@times,ddddJ_TdTdSigmadSigmadSigma,permute(dTdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdRdSigmadSigmadSigma,permute(dRdphi,[3,1,2])) ...
                            + bsxfun(@times,ddddJ_TdSigmadSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
                        ddddJ_TdphidphidTdT = bsxfun(@times,ddddJ_TdphidTdTdT,permute(dTdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidTdTdR,permute(dRdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidTdTdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                        ddddJ_TdphidphidTdR = bsxfun(@times,ddddJ_TdphidTdTdR,permute(dTdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidTdRdR,permute(dRdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidTdRdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                        ddddJ_TdphidphidRdR = bsxfun(@times,ddddJ_TdphidTdRdR,permute(dTdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidRdRdR,permute(dRdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidRdRdSigma,permute(dSigma_timedphi,[3,1,4,2]));
                        ddddJ_TdphidphidTdSigma = bsxfun(@times,ddddJ_TdphidTdTdSigma,permute(dTdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidTdRdSigma,permute(dRdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidTdSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                        ddddJ_TdphidphidRdSigma = bsxfun(@times,ddddJ_TdphidTdRdSigma,permute(dTdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidRdRdSigma,permute(dRdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidRdSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                        ddddJ_TdphidphidSigmadSigma = bsxfun(@times,ddddJ_TdphidTdSigmadSigma,permute(dTdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidRdSigmadSigma,permute(dRdphi,[3,1,4,2])) ...
                            + bsxfun(@times,ddddJ_TdphidSigmadSigmadSigma,permute(dSigma_timedphi,[3,1,4,2]));
                        ddddJ_TdphidphidphidT = bsxfun(@times,ddddJ_TdphidphidTdT,permute(dTdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidTdR,permute(dRdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidTdSigma,permute(dSigma_timedphi,[3,1,4,5,2]));
                        ddddJ_TdphidphidphidR = bsxfun(@times,ddddJ_TdphidphidTdR,permute(dTdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidRdR,permute(dRdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidRdSigma,permute(dSigma_timedphi,[3,1,4,5,2]));
                        ddddJ_TdphidphidphidSigma = bsxfun(@times,ddddJ_TdphidphidTdSigma,permute(dTdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidRdSigma,permute(dRdphi,[3,1,4,5,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidSigmadSigma,permute(dSigma_timedphi,[3,1,4,5,2]));
                        ddddJ_Tdphidphidphidphi = ddddJ_Tdphidphidphidphi + squeeze(sum(bsxfun(@times,ddddJ_TdphidphidphidT,permute(dTdphi,[3,1,4,5,6,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidphidR,permute(dRdphi,[3,1,4,5,6,2])) ...
                            + bsxfun(@times,ddddJ_TdphidphidphidSigma,permute(dSigma_timedphi,[3,1,4,5,6,2])),2));
                        
                        ddddJ_Tdbdphidphidphi = permute(sum(bsxfun(@times,ddddJ_Tdphidphidphidphi,permute(dphidb,[1,3,4,5,2])),1),[5,2,3,4,1]);
                        ddddJ_Tdbdbdphidphi = permute(sum(bsxfun(@times,ddddJ_Tdbdphidphidphi,permute(dphidb,[3,1,4,5,2])),2),[1,5,3,4,2]) ...
                            + permute(sum(bsxfun(@times,dddJ_Tdphidphidphi,permute(ddphidbdb,[1,4,5,2,3])),1),[4,5,2,3,1]);
                        ddddJ_Tdbdbdbdphi = permute(sum(bsxfun(@times,ddddJ_Tdbdbdphidphi,permute(dphidb,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                            + permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(ddphidbdb,[4,1,5,2,3])),2),[1,4,5,3,2]);
                        ddddJ_Tdbdbdbdb = permute(sum(bsxfun(@times,ddddJ_Tdbdbdbdphi,permute(dphidb,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(ddphidbdb,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        FIM.dbdb = ddddJ_Ddbdbdbdb + ddddJ_Tdbdbdbdb + squeeze(ddddJ_bdbdbdbdb);
                        
                        %% FIM.dbdbeta
                        
                        ddddJ_Ddbdbdbdbeta = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbdphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        ddddJ_Tdbdbdbdbeta = permute(sum(bsxfun(@times,ddddJ_Tdbdbdbdphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(ddphidbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        FIM.dbdbeta = ddddJ_Ddbdbdbdbeta + ddddJ_Tdbdbdbdbeta;
                        
                        %% FIM.dbetadbeta
                        
                        ddddJ_Ddbdbdbetadphi = permute(sum(bsxfun(@times,ddddJ_Ddbdbdphidphi,permute(dphidbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdphidphi,permute(ddphidbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
                        ddddJ_Ddbdbdbetadbeta = permute(sum(bsxfun(@times,ddddJ_Ddbdbdbetadphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Ddbdbdphi,permute(ddphidbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        ddddJ_Tdbdbdbetadphi = permute(sum(bsxfun(@times,ddddJ_Tdbdbdphidphi,permute(dphidbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
                            + permute(sum(bsxfun(@times,dddJ_Tdbdphidphi,permute(ddphidbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
                        ddddJ_Tdbdbdbetadbeta = permute(sum(bsxfun(@times,ddddJ_Tdbdbdbetadphi,permute(dphidbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
                            + permute(sum(bsxfun(@times,dddJ_Tdbdbdphi,permute(ddphidbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
                        
                        FIM.dbetadbeta = ddddJ_Ddbdbdbetadbeta + ddddJ_Tdbdbdbetadbeta;
                        
                        %% FIM.dbddelta
                        
                        FIM.dbddelta = ddddJ_bdbdbdbddelta;
                        
                        %% FIM.ddeltaddelta
                        
                        FIM.ddeltaddelta = ddddJ_bdbdbddeltaddelta;
                        
                        %% FIM.dbetaddelta
                        FIM.dbetaddelta = zeros(length(b),length(b),length(beta),length(b));
                        
                    end
                end
            end
        end
    end
    
    if nargout >=1
        varargout{1} = J.val;
    end
    if nargout >= 2
        varargout{2} = J.db;
    end
    if nargout >= 3
        varargout{3} = FIM.val;
    end
    if nargout >= 4
        varargout{4} = J;
        varargout{5} = FIM;
    end
    
end