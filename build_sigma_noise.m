% build_sigma_noise is an auxiliary function for logL_CE_w_grad_2 to
% construct the standard deviation matrix for observables and adapt it to
% the size of the measurements
%
% USAGE:
% ======
% [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi,ddddSigma_noisedphidphidphidphi] = build_sigma_noise(phi,Ym,s,Model,ind_y)
%
% INPUTS:
% =======
% phi ... mixed effect parameter
% Ym ... measurements
% Model ... model definition
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


function [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi,ddddSigma_noisedphidphidphidphi] = build_sigma_noise(phi,Ym,s,Model,ind_y)

nt = size(Ym,1);
ny = size(Ym,2);
np = length(phi);

sigma_noise = Model.exp{s}.sigma_noise(phi);

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
    dsigma_noisedphi = Model.exp{s}.dsigma_noisedphi(phi);
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
        ddsigma_noisedphidphi = Model.exp{s}.ddsigma_noisedphidphi(phi);
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
        dddsigma_noisedphidphidphi = Model.exp{s}.dddsigma_noisedphidphidphi(phi);
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
        ddddsigma_noisedphidphidphidphi = Model.exp{s}.ddddsigma_noisedphidphidphidphi(phi);
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
        I = ((1:np)-1)*numel(Ym);
        ind_1 = bsxfun(@plus,ind_y,I);
        dSigma_noisedphi = reshape(dSNdphi(ind_1(:)),[numel(ind_y),np]);
        if nargout >= 3 % second order derivatives        
            II = bsxfun(@plus,I,permute(I*np,[1,3,2]));
            ind_2 = bsxfun(@plus,ind_y,II);
            ddSigma_noisedphidphi = reshape(ddSNdphidphi(ind_2(:)),[numel(ind_y),np,np]);
            if nargout >= 4 % third order derivatives
                III = bsxfun(@plus,II,permute(I*np^2,[1,3,4,2]));
                ind_3 = bsxfun(@plus,ind_y,III);
                dddSigma_noisedphidphidphi = reshape(dddSNdphidphidphi(ind_3(:)),[numel(ind_y),np,np,np]);
                if nargout >= 5 % fourth order derivatives
                    IIII = bsxfun(@plus,III,permute(I*np^3,[1,3,4,5,2]));
                    ind_4 = bsxfun(@plus,ind_y,IIII);
                    ddddSigma_noisedphidphidphidphi = reshape(ddddSNdphidphidphidphi(ind_4(:)),[numel(ind_y),np,np,np,np]);
                end
            end
        end
    end
end