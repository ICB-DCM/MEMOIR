% build_sigma_time is an auxiliary function for logL_CE_w_grad_2 to
% construct the standard deviation matrix for events and adapt it to
% the size of the observed events
%
% USAGE:
% ======
% [Sigma_time,dSigma_timedphi,ddSigma_timedphidphi,dddSigma_timedphidphidphi,ddddSigma_timedphidphidphidphi] = build_sigma_time(phi,Tm,s,Model,ind_t)
%
% INPUTS:
% =======
% phi ... mixed effect parameter
% Tm ... observed events
% Model ... model definition
% s ... index of considered experiment
% ind_t ... indexing of events
%
% Outputs:
% ========
% Sigma_time ... standard deviation matrix for observables
% dSigma_timedphi ... gradient of standard deviation matrix for observables
% ...
%
% 2015/04/14 Fabian Froehlich


function [Sigma_time,dSigma_timedphi,ddSigma_timedphidphi,dddSigma_timedphidphidphi,ddddSigma_timedphidphidphidphi] = build_sigma_time(phi,Tm,s,Model,ind_t)
sigma_time = Model.exp{s}.sigma_time(phi);

np = length(phi);
nt = size(Tm,1);
nr = size(Tm,2);

if(size(sigma_time,1) == nt)
    if(size(sigma_time,2) == 1)
        Sigma_time = repmat(sigma_time,[1,nr]);
    elseif(size(sigma,2) == nr)
        Sigma_time = sigma_time;
    else
        error('Incompatible size of sigma_time parametrisation!')
    end
elseif(size(sigma_time,2) == nr)
    if(size(sigma_time,1) == 1)
        Sigma_time = repmat(sigma_time,[nt,1]);
    else
        error('Incompatible size of sigma_time parametrisation!')
    end
elseif(and(size(sigma_time,1)==1,size(sigma_time,2)==1))
    Sigma_time = repmat(sigma_time,size(Tm));
else
    error('Incompatible size of sigma_time parametrisation!')
end

if nargout >= 2 % first order derivatives
    dsigma_timedphi = Model.exp{s}.dsigma_timedphi(phi);
    dSigma_timedphi = zeros(length(ind_t),np);
    
    if(size(dsigma_timedphi,1) == nt)
        if(size(dsigma_timedphi,2) == 1)
            dSTdphi = repmat(dsigma_timedphi,[1,nr,1]);
        elseif(size(dsigma_timedphi,2) == nr)
            dSTdphi = dsigma_timedphi;
        end
    elseif(size(dsigma_timedphi,2) == nr)
        if(size(dsigma_timedphi,1) == 1)
            dSTdphi = repmat(dsigma_timedphi,[nt,1,1]);
        end
    elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
        dSTdphi = repmat(dsigma_timedphi,[size(Tm),1]);
    end
    
    if nargout >= 3 % second order derivatives
        ddsigma_timedphidphi = Model.exp{s}.ddsigma_timedphidphi(phi);
        ddSigma_timedphidphi = zeros(length(ind_t),np,np);
        if(size(dsigma_timedphi,1) == nt)
            if(size(dsigma_timedphi,2) == 1)
                ddSTdphidphi = repmat(ddsigma_timedphidphi,[1,nr,1,1]);
            elseif(size(dsigma_timedphi,2) == nr)
                ddSTdphidphi = ddsigma_timedphidphi ;
            end
        elseif(size(dsigma_timedphi,2) == nr)
            if(size(dsigma_timedphi,1) == 1)
                ddSTdphidphi = repmat(ddsigma_timedphidphi,[nt,1,1]);
            end
        elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
            ddSTdphidphi = repmat(ddsigma_timedphidphi,[size(Tm),1,1]);
        end
    end
    
    if nargout >= 4 % third order derivatives
        dddsigma_timedphidphidphi = Model.exp{s}.dddsigma_timedphidphidphi(phi);
        dddSigma_timedphidphidphi = zeros(length(ind_t),np,np,np);
        if(size(dsigma_timedphi,1) == nt)
            if(size(dsigma_timedphi,2) == 1)
                dddSTdphidphidphi = repmat(dddsigma_timedphidphidphi,[1,nr,1,1,1]);
            elseif(size(dsigma_timedphi,2) == nr)
                dddSTdphidphidphi = dddsigma_timedphidphidphi;
            end
        elseif(size(dsigma_timedphi,2) == nr)
            if(size(dsigma_timedphi,1) == 1)
                dddSTdphidphidphi = repmat(dddsigma_timedphidphidphi,[nt,1,1,1]);
            end
        elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
            dddSTdphidphidphi = repmat(dddsigma_timedphidphidphi,[size(Tm),1,1,1]);
        end
    end
    
    if nargout >= 5 % fourth order derivatives
        ddddsigma_timedphidphidphidphi = Model.exp{s}.ddddsigma_timedphidphidphidphi(phi);
        ddddSigma_timedphidphidphidphi = zeros(length(ind_t),np,np,np,np);
        if(size(dsigma_timedphi,1) == nt)
            if(size(dsigma_timedphi,2) == 1)
                ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[1,nr,1,1,1]);
            elseif(size(dsigma_timedphi,2) == nr)
                ddddSTdphidphidphidphi = ddddsigma_timedphidphidphidphi;
            end
        elseif(size(dsigma_timedphi,2) == nr)
            if(size(dsigma_timedphi,1) == 1)
                ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[nt,1,1,1]);
            end
        elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
            ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[size(Tm),1,1,1]);
        end
    end
    
    if nargout >= 2 % first order derivatives
        tmp = reshape(dSTdphi,[numel(Tm),np]);
        dSigma_timedphi = tmp(ind_t,:);
        if nargout >= 3 % second order derivatives
            tmp = reshape(ddSTdphidphi,[numel(Tm),np,np]);
            ddSigma_timedphidphi = tmp(ind_t,:,:);
            if nargout >= 4 % third order derivatives
                tmp = reshape(dddSTdphidphidphi,[numel(Tm),np,np,np]);
                dddSigma_timedphidphidphi = tmp(ind_t,:,:,:);
                if nargout >= 5 % fourth order derivatives
                    tmp = reshape(ddddSTdphidphidphidphi,[numel(Tm),np,np,np,np]);
                    ddddSigma_timedphidphidphidphi = tmp(ind_t,:,:,:,:);
                end
            end
        end
    end
    
end