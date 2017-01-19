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


function [Sigma_time] = build_sigma_time(phi,Tm,s,Model,ind_t,nderiv)

sigma_time = Model.sigma_time(phi);

np = length(phi);
nt = size(Tm,1);
nr = size(Tm,2);

if(size(sigma_time,1) == nt)
    if(size(sigma_time,2) == 1)
        Sigma_time.val = repmat(sigma_time,[1,nr]);
    elseif(size(sigma,2) == nr)
        Sigma_time.val = sigma_time;
    else
        error('Incompatible size of sigma_time parametrisation!')
    end
elseif(size(sigma_time,2) == nr)
    if(size(sigma_time,1) == 1)
        Sigma_time.val = repmat(sigma_time,[nt,1]);
    else
        error('Incompatible size of sigma_time parametrisation!')
    end
elseif(and(size(sigma_time,1)==1,size(sigma_time,2)==1))
    Sigma_time.val = repmat(sigma_time,size(Tm));
else
    error('Incompatible size of sigma_time parametrisation!')
end

if nderiv >= 1 % first order derivatives
    dsigma_timedphi = Model.dsigma_timedphi(phi);
    
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
    
    if nderiv >= 2 % second order derivatives
        ddsigma_timedphidphi = Model.ddsigma_timedphidphi(phi);
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
    
    if nderiv >= 3 % third order derivatives
        dddsigma_timedphidphidphi = Model.dddsigma_timedphidphidphi(phi);
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
    
%     if nderiv >= 4 % fourth order derivatives
%         ddddsigma_timedphidphidphidphi = Model.ddddsigma_timedphidphidphidphi(phi);
%         if(size(dsigma_timedphi,1) == nt)
%             if(size(dsigma_timedphi,2) == 1)
%                 ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[1,nr,1,1,1]);
%             elseif(size(dsigma_timedphi,2) == nr)
%                 ddddSTdphidphidphidphi = ddddsigma_timedphidphidphidphi;
%             end
%         elseif(size(dsigma_timedphi,2) == nr)
%             if(size(dsigma_timedphi,1) == 1)
%                 ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[nt,1,1,1]);
%             end
%         elseif(and(size(dsigma_timedphi,1)==1,size(dsigma_timedphi,2)==1))
%             ddddSTdphidphidphidphi = repmat(ddddsigma_timedphidphidphidphi,[size(Tm),1,1,1]);
%         end
%     end
%     
    if nderiv >= 1 % first order derivatives
        Sigma_time.dphi = reshape(dSTdphi(ind_t,:,:),[sum(ind_t)*size(Tm,2),np]);
        if nderiv >= 2 % second order derivatives    
            Sigma_time.dphidphi = reshape(ddSTdphidphi(ind_t,:,:,:),[sum(ind_t)*size(Tm,2),np,np]);
            if nderiv >= 3 % third order derivatives
                Sigma_time.dphidphidphi = reshape(dddSTdphidphidphi(ind_t,:,:,:,:),[sum(ind_t)*size(Tm,2),np,np,np]);
%                 if nargout >= 5 % fourth order derivatives
%                     ddddSigma_timedphidphidphidphi = reshape(ddddSTdphidphidphidphi(ind_t,:,:,:,:,:),[sum(ind_t)*size(Tm,2),np,np,np,np]);
%                 end
            end
        end
    end
    
end