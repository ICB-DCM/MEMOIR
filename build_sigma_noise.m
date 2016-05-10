% build_sigma_noise is an auxiliary function for logL_CE_w_grad_2 to
% construct the standard deviation matrix for observables and adapt it to
% the size of the measurements
%
% USAGE:
% ======
% [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi,ddddSigma_noisedphidphidphidphi] = build_sigma_noise(phi,Ym,s,model,ind_y)
%
% INPUTS:
% =======
% phi ... mixed effect parameter
% Ym ... measurements
% model ... model definition
% s ... index of considered experiment
% ind_y ... indexing of observables
%
% Outputs:
% ========
% Sigma_noise ... standard deviation matrix for observables
% dSigma_noisedphi ... gradient of standard deviation matrix for observables
% ...
%
% 2015/04/14 Fabian Froehlich


function [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi,ddddSigma_noisedphidphidphidphi] = build_sigma_noise(phi,Ym,s,model,ind_y)

nt = size(Ym,1);
ny = size(Ym,2);
np = length(phi);

sigma_noise = model.sigma_noise(phi);

if(size(sigma_noise,1) == nt)
    if(size(sigma_noise,2) == 1)
        Sigma_noise = repmat(sigma_noise,[1,ny]);
    elseif(size(sigma,2) == ny)
        Sigma_noise = sigma_noise;
    else
        error('Incompatible size of sigma_noise parametrisation!')
    end
elseif(size(sigma_noise,2) == ny)
    if(size(sigma_noise,1) == 1)
        Sigma_noise = repmat(sigma_noise,[nt,1]);
    else
        error('Incompatible size of sigma_noise parametrisation!')
    end
elseif(and(size(sigma_noise,1)==1,size(sigma_noise,2)==1))
    Sigma_noise = repmat(sigma_noise,size(Ym));
else
    error('Incompatible size of sigma_noise parametrisation!')
end

if nargout >= 2 % first order derivatives
    dsigma_noisedphi = model.dsigma_noisedphi(phi);
    dSigma_noisedphi = zeros(length(ind_y),np);
    
    if(size(dsigma_noisedphi,1) == nt)
        if(size(dsigma_noisedphi,2) == 1)
            dSNdphi = repmat(dsigma_noisedphi,[1,ny,1]);
        elseif(size(dsigma_noisedphi,2) == ny)
            dSNdphi = dsigma_noisedphi;
        end
    elseif(size(dsigma_noisedphi,2) == ny)
        if(size(dsigma_noisedphi,1) == 1)
            dSNdphi = repmat(dsigma_noisedphi,[nt,1,1]);
        end
    elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
        dSNdphi = repmat(dsigma_noisedphi,[size(Ym),1]);
    end
    
    if nargout >= 3 % second order derivatives
        ddsigma_noisedphidphi = model.ddsigma_noisedphidphi(phi);
        ddSigma_noisedphidphi = zeros(length(ind_y),np,np);
        if(size(dsigma_noisedphi,1) == nt)
            if(size(dsigma_noisedphi,2) == 1)
                ddSNdphidphi = repmat(ddsigma_noisedphidphi,[1,ny,1,1]);
            elseif(size(dsigma_noisedphi,2) == ny)
                ddSNdphidphi = ddsigma_noisedphidphi ;
            end
        elseif(size(dsigma_noisedphi,2) == ny)
            if(size(dsigma_noisedphi,1) == 1)
                ddSNdphidphi = repmat(ddsigma_noisedphidphi,[nt,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            ddSNdphidphi = repmat(ddsigma_noisedphidphi,[size(Ym),1,1]);
        end
    end
    
    if nargout >= 4 % third order derivatives
        dddsigma_noisedphidphidphi = model.dddsigma_noisedphidphidphi(phi);
        dddSigma_noisedphidphidphi = zeros(length(ind_y),np,np,np);
        if(size(dsigma_noisedphi,1) == nt)
            if(size(dsigma_noisedphi,2) == 1)
                dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[1,ny,1,1,1]);
            elseif(size(dsigma_noisedphi,2) == ny)
                dddSNdphidphidphi = dddsigma_noisedphidphidphi;
            end
        elseif(size(dsigma_noisedphi,2) == ny)
            if(size(dsigma_noisedphi,1) == 1)
                dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[nt,1,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[size(Ym),1,1,1]);
        end
    end
    
    if nargout >= 5 % fourth order derivatives
        ddddsigma_noisedphidphidphidphi = model.ddddsigma_noisedphidphidphidphi(phi);
        ddddSigma_noisedphidphidphidphi = zeros(length(ind_y),np,np,np,np);
        if(size(dsigma_noisedphi,1) == nt)
            if(size(dsigma_noisedphi,2) == 1)
                ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[1,ny,1,1,1]);
            elseif(size(dsigma_noisedphi,2) == ny)
                ddddSNdphidphidphidphi = ddddsigma_noisedphidphidphidphi;
            end
        elseif(size(dsigma_noisedphi,2) == ny)
            if(size(dsigma_noisedphi,1) == 1)
                ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[nt,1,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[size(Ym),1,1,1]);
        end
    end
    
    if nargout >= 2 % first order derivatives
        dSigma_noisedphi = reshape(dSNdphi(ind_y,:,:),[sum(ind_y)*size(Ym,2),np]);
        if nargout >= 3 % second order derivatives    
            ddSigma_noisedphidphi = reshape(ddSNdphidphi(ind_y,:,:,:),[sum(ind_y)*size(Ym,2),np,np]);
            if nargout >= 4 % third order derivatives
                dddSigma_noisedphidphidphi = reshape(dddSNdphidphidphi(ind_y,:,:,:,:),[sum(ind_y)*size(Ym,2),np,np,np]);
                if nargout >= 5 % fourth order derivatives
                    ddddSigma_noisedphidphidphidphi = reshape(ddddSNdphidphidphidphi(ind_y,:,:,:,:,:),[sum(ind_y)*size(Ym,2),np,np,np,np]);
                end
            end
        end
    end
end