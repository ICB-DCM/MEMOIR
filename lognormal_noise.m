%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_D = lognormal_noise(Y,Ym,Sigma,ind,nderiv)

    % linearize inputs
    Ym = Ym(ind,:);
    Sigma = Sigma(ind,:);
    
    if nderiv >=0
        % J_D
        J_D.val = nansum(0.5*((log(Y(:)) - Ym(:))./Sigma(:)).^2 + 0.5*log(sqrt(2*pi)*(Sigma(:).^2).*(Y(:).^2)));
        if nderiv >= 1
            % dJ_DdY
            J_D.dY = transpose((log(Y(:)) - Ym(:))./(Sigma(:).^2).*(1./Y) + 1./Y(:));
            % dJ_DdSigma
            J_D.dSigma = transpose(- (((log(Y(:)) - Ym(:)).^2)./(Sigma(:).^3)) + 1./Sigma(:));
            if nderiv >= 2
                %ddJ_DdYdY
                J_D.dYdY = (1-(log(Y(:)) - Ym(:)))./(Sigma(:).^2.*Y.^2) - 1./(Y(:).^2);
                %ddJ_DdYdSigma
                J_D.dYdSigma = -2*(Y(:) - Ym(:))./(Sigma(:).^3.*Y);
                %ddJ_DdSigmadSigma
                J_D.dSigmadSigma = 3*(((Y(:) - Ym(:)).^2)./(Sigma(:).^4)) - 1./(Y(:).^2);
                if nderiv >= 3
                    %dddJ_DdYdYdY
                    J_D.dYdYdY = - 3./(Sigma(:).^2.*Y.^3)-2./(Y.^3);
                    %dddJ_DdYdYdSigma
                    J_D.dYdYdSigma = (1-(log(Y(:)) - Ym(:)))./(Sigma(:).^3.*Y.^2);
                    %dddJ_DdYdSigmadSigma
                    J_D.dYdSigmadSigma = + 6*(Y(:) - Ym(:))./(Sigma(:).^4.*Y);
                    %dddJ_DdSigmadSigmadSigma
                    J_D.dSigmadSigmadSigma = - 12*(((Y(:) - Ym(:)).^2)./(Sigma(:).^5)) + 2./(Sigma(:).^3);
                    if nderiv >= 4
                        %ddddJ_DdYdYdYdY
                        J_D.dYdYdYdY = +9./(Sigma(:).^2.*Y.^4)+6./(Y.^4);
                        %ddddJ_DdYdYdYdSigma
                        J_D.dYdYdYdSigma = 6./(Sigma(:).^3.*Y.^3)-2./(Y.^3);
                        %ddddJ_DdYdYdSigmadSigma
                        J_D.dYdYdSigmadSigma = -3*(1-(log(Y(:)) - Ym(:)))./(Sigma(:).^4.*Y.^2);
                        %ddddJ_DdYdSigmadSigmadSigma
                        J_D.dYdSigmadSigmadSigma = - 24*(Y(:) - Ym(:))./(Sigma(:).^5.*Y);
                        %ddddJ_DdSigmadSigmadSigmadSigma
                        J_D.dSigmadSigmadSigmadSigma = - 60*(((Y(:) - Ym(:)).^2)./(Sigma(:).^6)) + 6./(Sigma(:).^4);
                    end
                end
            end
        end
    end
end