% simulate_trajectory is an auxiliary function for logL_CE_w_grad_2 to
% simulate individual trajectories for SCTL data
%
% USAGE:
% ======
% [Y,T,R,dYdphi,dTdphi,dRdphi,ddYdphidphi,ddTdphidphi,ddRdphidphi] = simulate_trajectory(t,phi,Model,kappa,s,ind_t,ind_y)
%
% INPUTS:
% =======
% t ... time vector
% phi ... mixed effect parameter
% Model ... model definition, here only the .model and .integration 
%     field is required
% kappa ... experimental condition
% s ... index of considered experiment
% ind_t ... indexing of events
% ind_y ... indexing of observables
%
% Outputs:
% ========
% Y ... model output of observables
% T ... model output of event (zero-crossings of rootfunction)
% R ... model output of rootvalue (value of rootfunction at time of event)
%     this should only be nonzero for T = max(t) i.e. events that did not
%     happen
% dYdphi ... sensitivity of observable
% dTdphi ... sensitivity of events
% dRdphi ... sensitivity of rootfunction at events
% ddYdphidphi ... second order sensitivity of observable
% ddTdphidphi ... second order sensitivity of events
% ddRdphidphi ... second order sensitivity of rootfunction at events
%
% 2015/04/14 Fabian Froehlich

function [ YY,TT,RR ] = simulate_trajectory(t,phi,Model,kappa,s,ind_t,ind_y,nderiv)

optionmu.atol = 1e-8;
optionmu.rtol = 1e-12;
optionmu.nmaxevent = length(ind_t)+10;
optionmu.maxsteps = 1e5;

if nderiv == 0
    optionmu.sensi = 0; % number of requested sensitivities
    sol = Model.model(t,phi,kappa,optionmu);
elseif nderiv == 1
    optionmu.sensi = 1; % number of requested sensitivities
    sol = Model.model(t,phi,kappa,optionmu);
    dYdphi = sol.sy(1:(end-(sum(ind_y)<length(t))),:,:);
    if(isfield(sol,'sz') && sum(ind_t)>0)
        dTdphi = sol.sz(ind_t,:,:);
        dRdphi = sol.srz(ind_t,:,:);
    else
        dTdphi = zeros(0,sum(ind_t),length(phi));
        dRdphi = zeros(0,sum(ind_t),length(phi));
    end
else
    optionmu.sensi = 2; % number of requested sensitivities
    optionmu.linsol = 9;
    sol = Model.model(t,phi,kappa,optionmu);
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(phi,@(phi) Model.model(t,phi,kappa,optionmu),1e-5,'y','sy')
    dYdphi = sol.sy(1:(end-(sum(ind_y)<length(t))),:,:);
    if(isfield(sol,'sz') && sum(ind_t)>0)
        dTdphi = sol.sz(ind_t,:,:);
        dRdphi = sol.srz(ind_t,:,:);
    else
        dTdphi = zeros(0,sum(ind_t),length(phi));
        dRdphi = zeros(0,sum(ind_t),length(phi));
    end
    try
        ddYdphidphi = sol.s2y(1:(end-(sum(ind_y)<length(t))),:,:,:);
        ddTdphidphi = sol.s2z(ind_t,:,:,:);
        ddRdphidphi = sol.s2rz(ind_t,:,:,:);
    catch
        % fallback if no second order sensitivities are available
        ddYdphidphi = zeros(length(t),size(sol.sy,2),length(phi),length(phi));
        if(isfield(sol,'sz') && sum(ind_t)>0)
            ddTdphidphi = zeros(sum(ind_t),size(sol.sz,2),length(phi),length(phi));
            ddRdphidphi = zeros(sum(ind_t),size(sol.sz,2),length(phi),length(phi));
        else
            ddTdphidphi = zeros(sum(ind_t),1,length(phi),length(phi));
            ddRdphidphi = zeros(sum(ind_t),1,length(phi),length(phi));
        end
        
    end
end
if(sol.status <0)
    error('integration failed');
end
Y = sol.y(1:(end-(sum(ind_y)<length(t))),:);
if(isfield(sol,'z'))
    TT.val = sol.z;
    RR.val = sol.rz;
else
    TT.val = zeros(0,sum(ind_t));
    RR.val = zeros(0,sum(ind_t));
end
YY.val = Y(:);

% Apply indexing

            
if nderiv >= 1
    YY.dphi = reshape(dYdphi,[size(dYdphi,1)*size(dYdphi,2),size(dYdphi,3)]);
    TT.dphi = reshape(dTdphi,[size(dTdphi,1)*size(dTdphi,2),size(dTdphi,3)]);
    RR.dphi = reshape(dRdphi,[size(dRdphi,1)*size(dRdphi,2),size(dRdphi,3)]);
    if nderiv >= 2
        YY.dphidphi = reshape(ddYdphidphi,[size(ddYdphidphi,1)*size(ddYdphidphi,2),size(ddYdphidphi,3),size(ddYdphidphi,4)]);
        TT.dphidphi = reshape(ddTdphidphi,[size(ddTdphidphi,1)*size(ddTdphidphi,2),size(ddTdphidphi,3),size(ddTdphidphi,4)]);
        RR.dphidphi = reshape(ddRdphidphi,[size(ddRdphidphi,1)*size(ddRdphidphi,2),size(ddRdphidphi,3),size(ddRdphidphi,4)]);
    end
end

end

