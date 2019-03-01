%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY(indnan)FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_D = normal_noise(Y,Ym,Sigma,ind,nderiv)
    
    % linearize inputs
    Ym = Ym(ind,:);
    Sigma = Sigma(ind,:);
    ind_Ym = repmat(ind, size(Ym, 2), 1);
    
    if nderiv >= 0
        % J_D
        J_D.val = nansum(0.5*((Y - Ym(:))./Sigma(:)).^2 + 0.5*log(2*pi*Sigma(:).^2));
        if nderiv >= 1
            % dJ_DdY
            J_D.dY = transpose((Y - Ym(:))./(Sigma(:).^2));
            J_D.dY(isnan(J_D.dY)) = 0;
            % dJ_DdSigma
            J_D.dSigma = transpose(- (((Y - Ym(:)).^2)./(Sigma(:).^3)) + 1./Sigma(:));
            J_D.dSigma(isnan(J_D.dSigma)) = 0;
            if nderiv >= 2
                %ddJ_DdYdY
                J_D.dYdY = diag(1./(Sigma(:).^2));
                J_D.dYdY(isnan(J_D.dYdY)) = 0;
                %ddJ_DdYdSigma
                J_D.dYdSigma = diag(-2*(Y - Ym(:))./(Sigma(:).^3));
                J_D.dYdSigma(isnan(J_D.dYdSigma)) = 0;
                %ddJ_DdSigmadSigma
                J_D.dSigmadSigma = diag(3*(((Y - Ym(:)).^2)./(Sigma(:).^4)) - 1./(Sigma(:).^2));
                J_D.dSigmadSigma(isnan(J_D.dSigmadSigma)) = 0;
                if nderiv >= 3
                    %dddJ_DdYdYdY
                    J_D.dYdYdY = zeros([sum(ind_Ym),sum(ind_Ym),sum(ind_Ym)]);
                    J_D.dYdYdY(isnan(J_D.dYdYdY)) = 0;
                    %dddJ_DdYdYdSigma
                    J_D.dYdYdSigma = bsxfun(@times,diag(- 2./(Sigma(:).^3)),permute(eye(sum(ind_Ym)),[3,1,2]));
                    J_D.dYdYdSigma(isnan(J_D.dYdYdSigma)) = 0;
                    %dddJ_DdYdSigmadSigma
                    J_D.dYdSigmadSigma = bsxfun(@times,diag(6*(Y - Ym(:))./(Sigma(:).^4)),permute(eye(sum(ind_Ym)),[3,1,2]));
                    J_D.dYdSigmadSigma(isnan(J_D.dYdSigmadSigma)) = 0;
                    %dddJ_DdSigmadSigmadSigma
                    J_D.dSigmadSigmadSigma = bsxfun(@times,diag(- 12*(((Y - Ym(:)).^2)./(Sigma(:).^5)) + 2./(Sigma(:).^3)),permute(eye(sum(ind_Ym)),[3,1,2]));
                    J_D.dSigmadSigmadSigma(isnan(J_D.dSigmadSigmadSigma)) = 0;
    %                 if nargout >= 11
    %                     %ddddJ_DdYdYdYdY
    %                     varargout{11} = transpose(zeros(size(Y)));
    %                     %ddddJ_DdYdYdYdSigma
    %                     varargout{12} = transpose(zeros(size(Y)));
    %                     %ddddJ_DdYdYdSigmadSigma
    %                     varargout{13} = transpose(6./(Sigma(:).^4));
    %                     %ddddJ_DdYdSigmadSigmadSigma
    %                     varargout{14} = transpose(- 24*((Y - Ym(:))./(Sigma(:).^5)));
    %                     %ddddJ_DdSigmadSigmadSigmadSigma
    %                     varargout{15} = transpose(60*(((Y - Ym(:)).^2)./(Sigma(:).^6)) - 6./(Sigma(:).^4));
    %                 end
                end
            end
        end
    end
end