function [J_D,J_T,Sim] = objective_phi(model,data,phi,s,i,options,nderiv,FIMflag)

Ym = data.SCTL.Y(:,:,i);
ind_y = data.SCTL.ind_y(:,i);
Tm  = data.SCTL.T(:,:,i);
ind_t = data.SCTL.ind_t(:,i);

% build standard deviation matrices
Sigma_noise = build_sigma_noise(phi,Ym,s,model,ind_y,nderiv);
Sigma_time = build_sigma_time(phi,Tm,s,model,ind_t,nderiv);

t = data.SCTL.time(ind_y);
if(t(end)<max(max(Tm))*1.2)
    t = [t;max(max(Tm))*1.2];
end

% simulate trajectory
try
    [Y,T,R] = simulate_trajectory(t,phi,model,data.condition,s,ind_t,ind_y,nderiv);
catch err
    % if simulation fails, return Inf
    J_D.val = Inf;
    J_T.val = Inf;
    if(nargout>=3)
        Sim.SCTL_Y = Y.val;
        Sim.SCTL_T = T.val(ind_t);
        Sim.SCTL_R = R.val(ind_t);
        Sim.SCTL_Sigma_Y = Sigma_noise.val(ind_y);
        Sim.SCTL_Sigma_T = Sigma_time.val(ind_t);
    end
    return
end

% noise model
if(~isfield(model,'noise_model'))
    model.noise_model = 'normal';
end

switch(model.noise_model)
    case 'normal'
        noisedist = @normal_noise_optims;
    case 'lognormal'
        noisedist = @lognormal_noise;
    case 'tdist'
        noisedist = @tdist_noise;
end
J_D = noisedist(Y.val,Ym,Sigma_noise.val,ind_y,nderiv+FIMflag);

% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Y.val,@(Y) noisedist(Y,Ym,Sigma_noise.val,ind_y,nderiv+(nderiv==1)),1e-6,'val','dY')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Y.val,@(Y) noisedist(Y,Ym,Sigma_noise.val,ind_y,nderiv+(nderiv==1)),1e-6,'dY','dYdY')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Y.val,@(Y) noisedist(Y,Ym,Sigma_noise.val,ind_y,nderiv+(nderiv==1)),1e-6,'dYdY','dYdYdY')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Sigma_noise.val,@(Sigma) noisedist(Y.val,Ym,Sigma,ind_y,nderiv+(nderiv==1)),1e-6,'val','dSigma')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Sigma_noise.val,@(Sigma) noisedist(Y.val,Ym,Sigma,ind_y,nderiv+(nderiv==1)),1e-6,'dSigma','dSigmadSigma')


% event model
if(~isfield(model,'time_model'))
    model.time_model = 'normal';
end

switch(model.time_model)
    case 'normal'
        timedist = @normal_time;
end
J_T = timedist(T.val,Tm,R.val,Sigma_time.val,ind_t,nderiv+FIMflag);
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(T.val,@(T) timedist(T,Tm,R.val,Sigma_time.val,ind_t,nderiv+(nderiv==1)),1e-6,'val','dT')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(T.val,@(T) timedist(T,Tm,R.val,Sigma_time.val,ind_t,nderiv+(nderiv==1)),1e-6,'dT','dTdT')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(T.val,@(T) timedist(T,Tm,R.val,Sigma_time.val,ind_t,nderiv+(nderiv==1)),1e-6,'dTdT','dTdTdT')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(R.val,@(R) timedist(T.val,Tm,R,Sigma_time.val,ind_t,nderiv+(nderiv==1)),1e-6,'val','dR')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(R.val,@(R) timedist(T.val,Tm,R,Sigma_time.val,ind_t,nderiv+(nderiv==1)),1e-6,'dR','dRdR')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(R.val,@(R) timedist(T.val,Tm,R,Sigma_time.val,ind_t,nderiv+(nderiv==1)),1e-6,'dRdR','dRdRdR')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Sigma_time.val,@(Sigma) timedist(T.val,Tm,R.val,Sigma,ind_t,nderiv+(nderiv==1)),1e-6,'val','dSigma')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Sigma_time.val,@(Sigma) timedist(T.val,Tm,R.val,Sigma,ind_t,nderiv+(nderiv==1)),1e-6,'dSigma','dSigmadSigma')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(Sigma_time.val,@(Sigma) timedist(T.val,Tm,R.val,Sigma,ind_t,nderiv+(nderiv==1)),1e-6,'dSigmadSigma','dSigmadSigmadSigma')



if nderiv >= 1
    J_D.dphi = chainrule(J_D.dY,Y.dphi) + chainrule(J_D.dSigma,Sigma_noise.dphi);
    
    J_T.dphi = chainrule(J_T.dT,T.dphi) + chainrule(J_T.dR,R.dphi) + chainrule(J_T.dSigma,Sigma_time.dphi);
    
    if(FIMflag || nderiv >=2 )
        J_D.dphidY = permute(sum(bsxfun(@times,J_D.dYdY,permute(Y.dphi,[1,3,2])) ...
            + bsxfun(@times,J_D.dYdSigma,permute(Sigma_noise.dphi,[3,1,2])),1),[3,2,1]);
        J_D.dphidSigma = permute(sum(bsxfun(@times,J_D.dYdSigma,permute(Y.dphi,[1,3,2])) ...
            + bsxfun(@times,J_D.dSigmadSigma,permute(Sigma_noise.dphi,[1,3,2])),1),[3,2,1]);
        
        J_D.FIM = permute(sum(bsxfun(@times,J_D.dphidY,permute(Y.dphi,[3,1,2])) ...
            + bsxfun(@times,J_D.dphidSigma,permute(Sigma_noise.dphi,[3,1,2])),2),[1,3,2]);
        
        J_T.dphidT = permute(sum(bsxfun(@times,J_T.dTdT,permute(T.dphi,[1,3,2])) ...
            + bsxfun(@times,J_T.dTdR,permute(R.dphi,[1,3,2])) ...
            + bsxfun(@times,J_T.dTdSigma,permute(Sigma_time.dphi,[3,1,2])),1),[3,2,1]);
        J_T.dphidR = permute(sum(bsxfun(@times,J_T.dTdR,permute(T.dphi,[1,3,2])) ...
            + bsxfun(@times,J_T.dRdR,permute(R.dphi,[1,3,2])) ...
            + bsxfun(@times,J_T.dRdSigma,permute(Sigma_time.dphi,[1,3,2])),1),[3,2,1]);
        J_T.dphidSigma = permute(sum(bsxfun(@times,J_T.dTdSigma,permute(T.dphi,[1,3,2])) ...
            + bsxfun(@times,J_T.dRdSigma,permute(R.dphi,[1,3,2])) ...
            + bsxfun(@times,J_T.dSigmadSigma,permute(Sigma_time.dphi,[1,3,2])),1),[3,2,1]);
        
        J_T.FIM = permute(sum(bsxfun(@times,J_T.dphidT,permute(T.dphi,[3,1,2])) ...
            + bsxfun(@times,J_T.dphidR,permute(R.dphi,[3,1,2])) ...
            + bsxfun(@times,J_T.dphidSigma,permute(Sigma_time.dphi,[3,1,2])),2),[1,3,2]);
    end
    if nderiv >= 2
        % exact
        J_D.dphidphi = J_D.FIM ...
            + permute(sum(bsxfun(@times,J_D.dY,permute(Y.dphidphi,[4,1,2,3])),2),[3,4,1,2]);
        
        J_T.dphidphi = J_T.FIM...
            + permute(sum(bsxfun(@times,J_T.dT,permute(T.dphidphi,[4,1,2,3])),2),[3,4,1,2]);
        
        if nderiv >= 3
            temp = permute(sum(bsxfun(@times,J_D.dphidY,permute(Y.dphidphi,[4,1,2,3])) ...
                + bsxfun(@times,J_D.dphidSigma,permute(Sigma_noise.dphidphi,[4,1,2,3])),2),[1,3,4,2]);
            
            J_D.FIMdphi = permute(chainrule(J_D.dSigma,Sigma_noise.dphidphidphi),[2,3,4,1]) ...
                + permute(temp,[2,1,3]) + permute(temp,[1,2,3]);
            
            J_D.dphidYdY = permute(sum(bsxfun(@times,J_D.dYdYdY,permute(Y.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_D.dYdYdSigma,permute(Sigma_noise.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_D.dphidYdSigma = permute(sum(bsxfun(@times,J_D.dYdYdSigma,permute(Y.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_D.dYdSigmadSigma,permute(Sigma_noise.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_D.dphidSigmadSigma = permute(sum(bsxfun(@times,J_D.dYdSigmadSigma,permute(Y.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_D.dSigmadSigmadSigma,permute(Sigma_noise.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_D.dphidphidY = permute(sum(bsxfun(@times,J_D.dphidYdY,permute(Y.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_D.dphidYdSigma,permute(Sigma_noise.dphi,[3,1,4,2])),2),[1,4,3,2]);
            J_D.dphidphidSigma = permute(sum(bsxfun(@times,J_D.dphidYdSigma,permute(Y.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_D.dphidSigmadSigma,permute(Sigma_noise.dphi,[3,1,4,2])),2),[1,4,3,2]);
            
            J_D.FIMdphi = J_D.FIMdphi + permute(sum(bsxfun(@times,J_D.dphidphidY,permute(Y.dphi,[3,4,1,2])) ...
                + bsxfun(@times,J_D.dphidphidSigma,permute(Sigma_noise.dphi,[3,4,1,2])),3),[1,2,4,3]);
            
            temp = permute(sum(bsxfun(@times,J_T.dphidT,permute(T.dphidphi,[4,1,2,3])) ...
                + bsxfun(@times,J_T.dphidR,permute(R.dphidphi,[4,1,2,3])) ...
                + bsxfun(@times,J_T.dphidSigma,permute(Sigma_time.dphidphi,[4,1,2,3])),2),[1,3,4,2]);
            
            J_T.FIMdphi = permute(chainrule(J_T.dSigma,Sigma_time.dphidphidphi),[2,3,4,1]) ...
                + permute(temp,[2,1,3]) + permute(temp,[1,2,3]);
            
            J_T.dphidTdT = permute(sum(bsxfun(@times,J_T.dTdTdT,permute(T.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dTdTdR,permute(R.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dTdTdSigma,permute(Sigma_time.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_T.dphidTdR = permute(sum(bsxfun(@times,J_T.dTdTdR,permute(T.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dTdRdR,permute(R.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dTdRdSigma,permute(Sigma_time.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_T.dphidRdR = permute(sum(bsxfun(@times,J_T.dTdRdR,permute(T.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dRdRdR,permute(R.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dRdRdSigma,permute(Sigma_time.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_T.dphidTdSigma = permute(sum(bsxfun(@times,J_T.dTdTdSigma,permute(T.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dTdRdSigma,permute(R.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dTdSigmadSigma,permute(Sigma_time.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_T.dphidRdSigma = permute(sum(bsxfun(@times,J_T.dTdRdSigma,permute(T.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dRdRdSigma,permute(R.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dRdSigmadSigma,permute(Sigma_time.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_T.dphidSigmadSigma = permute(sum(bsxfun(@times,J_T.dTdSigmadSigma,permute(T.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dRdSigmadSigma,permute(R.dphi,[1,3,4,2])) ...
                + bsxfun(@times,J_T.dSigmadSigmadSigma,permute(Sigma_time.dphi,[1,3,4,2])),1),[4,2,3,1]);
            J_T.dphidphidT = permute(sum(bsxfun(@times,J_T.dphidTdT,permute(T.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_T.dphidTdR,permute(R.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_T.dphidTdSigma,permute(Sigma_time.dphi,[3,1,4,2])),2),[1,4,3,2]);
            J_T.dphidphidR = permute(sum(bsxfun(@times,J_T.dphidTdR,permute(T.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_T.dphidRdR,permute(R.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_T.dphidRdSigma,permute(Sigma_time.dphi,[3,1,4,2])),2),[1,4,3,2]);
            J_T.dphidphidSigma = permute(sum(bsxfun(@times,J_T.dphidTdSigma,permute(T.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_T.dphidRdSigma,permute(R.dphi,[3,1,4,2])) ...
                + bsxfun(@times,J_T.dphidSigmadSigma,permute(Sigma_time.dphi,[3,1,4,2])),2),[1,4,3,2]);
            
            J_T.FIMdphi = J_T.FIMdphi + permute(sum(bsxfun(@times,J_T.dphidphidT,permute(T.dphi,[3,4,1,2])) ...
                + bsxfun(@times,J_T.dphidphidR,permute(R.dphi,[3,4,1,2])) ...
                + bsxfun(@times,J_T.dphidphidSigma,permute(Sigma_time.dphi,[3,4,1,2])),3),[1,2,4,3]);
            
        end
    end
end

if(nargout>=3)
    Sim.SCTL_Y = Y.val;
    Sim.SCTL_T = T.val(ind_t);
    Sim.SCTL_R = R.val(ind_t);
    Sim.SCTL_Sigma_Y = Sigma_noise.val(ind_y);
    Sim.SCTL_Sigma_T = Sigma_time.val(ind_t);
end

end



