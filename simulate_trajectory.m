% simulate_trajectory is an auxiliary function for logL_CE_w_grad_2 to
% simulate individual trajectories for SCTL data
%
% USAGE:
% ======
% [Y,T,R,dYdphi,dTdphi,dRdphi,ddYdphidphi,ddTdphidphi,ddRdphidphi] = simulate_trajectory(t,phi,Model,kappa,s,ind_t)
%
% INPUTS:
% =======
% t ... time vector
% phi ... mixed effect parameter
% Model ... model definition, here only the .exp{s}.model and .integration 
%     field is required
% kappa ... experimental condition
% s ... index of considered experiment
% ind_t ... indexing of event
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

function [ Y,T,R, dYdphi,dTdphi,dRdphi, ddYdphidphi,ddTdphidphi,ddRdphidphi ] = simulate_trajectory(t,phi,Model,kappa,s,ind_t)

if nargout < 4
    option_simu.sensi = 0; % number of requested sensitivities
    sol = Model.exp{s}.model(t,phi,kappa,option_simu);
elseif nargout < 7
    option_simu.sensi = 1; % number of requested sensitivities
    sol = Model.exp{s}.model(t,phi,kappa,option_simu);
    dYdphi = sol.sy;
    dTdphi = sol.sroot(ind_t,:,:);
    dRdphi = sol.srootval(ind_t,:,:);
else
    option_simu.sensi = 2; % number of requested sensitivities
    sol = Model.exp{s}.model(t,phi,kappa,option_simu);
    dYdphi = sol.sy;
    dTdphi = sol.sroot(ind_t,:,:);
    dRdphi = sol.srootval(ind_t,:,:);
    try
        ddYdphidphi = sol.s2y;
        ddTdphidphi = sol.s2root(ind_t,:,:,:);
        ddRdphidphi = sol.s2rootval(ind_t,:,:,:);
    catch
        % fallback if no second order sensitivities are available
        ddYdphidphi = zeros(length(t),size(sol.sy,2),length(phi),length(phi));
        ddTdphidphi = zeros(sum(ind_t),size(sol.sroot,2),length(phi),length(phi));
        ddRdphidphi = zeros(sum(ind_t),size(sol.srootval,2),length(phi),length(phi));
    end
end
Y = sol.y;
T = sol.root;
R = sol.rootval;
end

