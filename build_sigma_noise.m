function [Sigma_noise,dSigma_noisedphi,ddSigma_noisedphidphi,dddSigma_noisedphidphidphi,ddddSigma_noisedphidphidphidphi] = build_sigma_noise(phi,Ym,s,Model,ind_y)

sigma_noise = Model.exp{s}.sigma_noise(phi);

if(size(sigma_noise,1) == size(Ym,1))
    if(size(sigma_noise,2) == 1)
        Sigma_noise = repmat(sigma_noise,[1,size(Ym,2)]);
    elseif(size(sigma,2) == size(Ym,2))
        Sigma_noise = sigma_noise;
    else
        error('Incompatible size of sigma_noise parametrisation!')
    end
elseif(size(sigma_noise,2) == size(Ym,2))
    if(size(sigma_noise,1) == 1)
        Sigma_noise = repmat(sigma_noise,[size(Ym,1),1]);
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
    dSigma_noisedphi = zeros(length(ind_y),length(phi));
    
    if(size(dsigma_noisedphi,1) == size(Ym,1))
        if(size(dsigma_noisedphi,2) == 1)
            dSNdphi = repmat(dsigma_noisedphi,[1,size(Ym,2),1]);
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            dSNdphi = dsigma_noisedphi;
        end
    elseif(size(dsigma_noisedphi,2) == size(Ym,2))
        if(size(dsigma_noisedphi,1) == 1)
            dSNdphi = repmat(dsigma_noisedphi,[size(Ym,1),1,1]);
        end
    elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
        dSNdphi = repmat(dsigma_noisedphi,[size(Ym),1]);
    end
    
    if nargout >= 3 % second order derivatives
        ddsigma_noisedphidphi = Model.exp{s}.ddsigma_noisedphidphi(phi);
        ddSigma_noisedphidphi = zeros(length(ind_y),length(phi),length(phi));
        if(size(dsigma_noisedphi,1) == size(Ym,1))
            if(size(dsigma_noisedphi,2) == 1)
                ddSNdphidphi = repmat(ddsigma_noisedphidphi,[1,size(Ym,2),1,1]);
            elseif(size(dsigma_noisedphi,2) == size(Ym,2))
                ddSNdphidphi = ddsigma_noisedphidphi ;
            end
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            if(size(dsigma_noisedphi,1) == 1)
                ddSNdphidphi = repmat(ddsigma_noisedphidphi,[size(Ym,1),1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            ddSNdphidphi = repmat(ddsigma_noisedphidphi,[size(Ym),1,1]);
        end
    end
    
    if nargout >= 4 % third order derivatives
        dddsigma_noisedphidphidphi = Model.exp{s}.dddsigma_noisedphidphidphi(phi);
        dddSigma_noisedphidphidphi = zeros(length(ind_y),length(phi),length(phi),length(phi));
        if(size(dsigma_noisedphi,1) == size(Ym,1))
            if(size(dsigma_noisedphi,2) == 1)
                dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[1,size(Ym,2),1,1,1]);
            elseif(size(dsigma_noisedphi,2) == size(Ym,2))
                dddSNdphidphidphi = dddsigma_noisedphidphidphi;
            end
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            if(size(dsigma_noisedphi,1) == 1)
                dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[size(Ym,1),1,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            dddSNdphidphidphi = repmat(dddsigma_noisedphidphidphi,[size(Ym),1,1,1]);
        end
    end
    
    if nargout >= 5 % fourth order derivatives
        ddddsigma_noisedphidphidphidphi = Model.exp{s}.ddddsigma_noisedphidphidphidphi(phi);
        ddddSigma_noisedphidphidphidphi = zeros(length(ind_y),length(phi),length(phi),length(phi),length(phi));
        if(size(dsigma_noisedphi,1) == size(Ym,1))
            if(size(dsigma_noisedphi,2) == 1)
                ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[1,size(Ym,2),1,1,1]);
            elseif(size(dsigma_noisedphi,2) == size(Ym,2))
                ddddSNdphidphidphidphi = ddddsigma_noisedphidphidphidphi;
            end
        elseif(size(dsigma_noisedphi,2) == size(Ym,2))
            if(size(dsigma_noisedphi,1) == 1)
                ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[size(Ym,1),1,1,1]);
            end
        elseif(and(size(dsigma_noisedphi,1)==1,size(dsigma_noisedphi,2)==1))
            ddddSNdphidphidphidphi = repmat(ddddsigma_noisedphidphidphidphi,[size(Ym),1,1,1]);
        end
    end
    
    for k = 1:length(phi) % first order derivatives
        temp = dSNdphi(:,:,k);
        dSigma_noisedphi(:,k) = temp(ind_y);
        if nargout >= 3 % second order derivatives
            for l = 1:length(phi)
                temp = ddSNdphidphi(:,:,k,l);
                ddSigma_noisedphidphi(:,k,l) = temp(ind_y);
                if nargout >= 4 % third order derivatives
                    for m = 1:length(phi)
                        temp = dddSNdphidphidphi(:,:,k,l,m);
                        dddSigma_noisedphidphidphi(:,k,l,m) = temp(ind_y);
                        if nargout >= 5 % fourth order derivatives
                            for n = 1:length(phi)
                                temp = ddddSNdphidphidphidphi(:,:,k,l,m,n);
                                ddddSigma_noisedphidphidphidphi(:,k,l,m,n) = temp(ind_y);
                            end
                        end
                    end
                end
            end
        end
    end
end