%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY(indnan)FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_D = normal_noise(Y,Ym,Sigma,ind,nderiv)

Ym = Ym(ind,:);
Sigma = Sigma(ind,:);

indnan = ~isnan(Ym(:));
res = (Y(indnan)- Ym(indnan))./Sigma(indnan);
    
if nderiv >= 0
    % J_D
    J_D.val = nansum(0.5*(res.^2 + log(2*pi*Sigma(indnan).^2)));
    if nderiv >= 1
        % dJ_DdY
        J_D.dY = zeros(1,length(Y));
        J_D.dY(indnan) = transpose(res./Sigma(indnan));
        % dJ_DdSigma
        J_D.dSigma = zeros(1,length(Y));
        J_D.dSigma(indnan) = transpose((1-res.^2)./Sigma(indnan));
        if nderiv >= 2
            %ddJ_DdYdY
            J_D.dYdY = zeros(length(Y));
            J_D.dYdY(indnan,indnan) = diag(1./(Sigma(indnan).^2));
            %ddJ_DdYdSigma
            J_D.dYdSigma = zeros(length(Y));
            J_D.dYdSigma(indnan,indnan) = diag(-2*res./(Sigma(indnan).^2));
            %ddJ_DdSigmadSigma
            J_D.dSigmadSigma = zeros(length(Y));
            J_D.dSigmadSigma(indnan,indnan) = diag((3*res.^2-1)./(Sigma(indnan).^2));
            if nderiv >= 3
                %dddJ_DdYdYdY
                J_D.dYdYdY = zeros([length(Y),length(Y),length(Y)]);
                %dddJ_DdYdYdSigma
                J_D.dYdYdSigma = zeros([length(Y),length(Y),length(Y)]);
                J_D.dYdYdSigma(indnan,indnan,indnan) = bsxfun(@times,diag(- 2./(Sigma(indnan).^3)),permute(eye(sum(indnan)),[3,1,2]));
                %dddJ_DdYdSigmadSigma
                J_D.dYdSigmadSigma = zeros([length(Y),length(Y),length(Y)]);
                J_D.dYdSigmadSigma(indnan,indnan,indnan) = bsxfun(@times,diag(6*res./(Sigma(indnan).^3)),permute(eye(sum(indnan)),[3,1,2]));
                %dddJ_DdSigmadSigmadSigma
                J_D.dSigmadSigmadSigma = zeros([length(Y),length(Y),length(Y)]);
                J_D.dSigmadSigmadSigma(indnan,indnan,indnan) = bsxfun(@times,diag((2+12*res.^2)./(Sigma(indnan).^3)),permute(eye(sum(indnan)),[3,1,2]));
%                 if nargout >= 11
%                     %ddddJ_DdYdYdYdY
%                     varargout{11} = transpose(zeros(size(Y)));
%                     %ddddJ_DdYdYdYdSigma
%                     varargout{12} = transpose(zeros(size(Y)));
%                     %ddddJ_DdYdYdSigmadSigma
%                     varargout{13} = transpose(6./(Sigma(indnan).^4));
%                     %ddddJ_DdYdSigmadSigmadSigma
%                     varargout{14} = transpose(- 24*(res./(Sigma(indnan).^5)));
%                     %ddddJ_DdSigmadSigmadSigmadSigma
%                     varargout{15} = transpose(60*((res.^2)./(Sigma(indnan).^6)) - 6./(Sigma(indnan).^4));
%                 end
            end
        end
    end
end
end