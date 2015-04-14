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
    
    if nargout >= 4 % third order derivatives
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
    
    if nargout >= 5 % fourth order derivatives
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
                if nargout >= 4 % third order derivatives
                    for m = 1:length(phi)
                        temp = dddSTdphidphidphi(:,:,k,l,m);
                        dddSigma_timedphidphidphi(:,k,l,m) = temp(ind_t);
                        if nargout >= 5 % fourth order derivatives
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