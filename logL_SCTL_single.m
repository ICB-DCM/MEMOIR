function [ varargout ] = logL_SCTL_single(phi,model,data,s,i)
%LOGL_SCTL_ Summary of this function goes here
%   Detailed explanation goes here

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
        
        ddJ_Dappdbdphi = transpose(squeeze(sum(bsxfun(@times,ddJ_Dappdphidphi,permute(dphidb,[3,1,2])),2)));
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
        
        ddJ_Tappdbdphi = transpose(squeeze(sum(bsxfun(@times,ddJ_Tappdphidphi,permute(dphidb,[3,1,2])),2)));
        ddJ_Tappdbdb = squeeze(sum(bsxfun(@times,ddJ_Tappdbdphi,permute(dphidb,[3,1,2,4])),2)) + chainrule(dJ_Tdphi,ddphidbdb);
        
        % FIM
        FIM.val = squeeze(ddJ_Dappdbdb) + squeeze(ddJ_Tappdbdb);
        
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
    
    
end
