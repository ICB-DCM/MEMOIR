%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_D = normal_noise(Y,Ym,Sigma,ind,nderiv)
Ym = Ym(ind,:);
Sigma = Sigma(ind,:);

S2 = nansum(((Y - Ym(:))).^2);

if nderiv >= 0
    % J_D
    J_D.val = 0.5*sum(ind)*(log(2*pi*S2/sum(ind)) + sum(ind));
    if nderiv >= 1
        % dJ_DdY
        J_D.dY = transpose(sum(ind)*(Y-Ym(:))/S2);
        % dJ_DdSigma
        J_D.dSigma = zeros(1,sum(ind));
        if nderiv >= 2
            %ddJ_DdYdY
            J_D.dYdY = sum(ind)*(diag(repmat(1/S2,[sum(ind),1]))-2*(Y-Ym(:))*transpose((Y-Ym(:)))/(S2^2));
            %ddJ_DdYdSigma
            J_D.dYdSigma = zeros([sum(ind),sum(ind)]);
            %ddJ_DdSigmadSigma
            J_D.dSigmadSigma = zeros([sum(ind),sum(ind)]);
            if nderiv >= 3
                %dddJ_DdYdYdY
                tmp = -2*bsxfun(@times,diag(repmat(1/(S2^2),[sum(ind),1])),permute((Y-Ym(:)),[3,2,1]));
                J_D.dYdYdY = sum(ind)*(bsxfun(@times,8*(Y-Ym(:))*transpose((Y-Ym(:))),permute((Y-Ym(:)),[3,2,1]))/(S2^3) + tmp + permute(tmp,[1,3,2]) + permute(tmp,[3,1,2]) - 0*bsxfun(@times,diag((Y-Ym(:))/(S2^2)),permute(eye(sum(ind)),[3,2,1])));
                %dddJ_DdYdYdSigma
                J_D.dYdYdSigma = zeros([sum(ind),sum(ind),sum(ind)]);
                %dddJ_DdYdSigmadSigma
                J_D.dYdSigmadSigma = zeros([sum(ind),sum(ind),sum(ind)]);
                %dddJ_DdSigmadSigmadSigma
                J_D.dSigmadSigmadSigma = zeros([sum(ind),sum(ind),sum(ind)]);
            end
        end
    end
end
end