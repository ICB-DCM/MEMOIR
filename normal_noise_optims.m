%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_D = normal_noise_optims(Y,Ym,Sigma,ind,nderiv)
Ym = Ym(ind,:);
Sigma = Sigma(ind,:);

N_i = sum(~isnan(Ym));
N = repmat(N_i,[size(Ym,1),1]);
S2_i = nansum((reshape(Y,size(Ym)) - Ym).^2)./sum(~isnan(Ym));
S2 = repmat(S2_i,[size(Ym,1),1]);


res = (Y-Ym(:))./S2(:);
res(isnan(res)) = 0;



if nderiv >= 0
    % J_D
    J_D.val = 0.5*nansum(~isnan(Ym(:)).*(log(2*pi*S2(:))+1));
    if nderiv >= 1
        % dJ_DdY
        J_D.dY = transpose(res);
        % dJ_DdSigma
        J_D.dSigma = zeros(1,length(Y));
        if nderiv >= 2
            %ddJ_DdYdY
            S2(isnan(Ym)) = Inf;
            resprod = (res)*transpose(res);
            krontmp = kron(diag(1./N_i),ones(size(Ym,1)));
            J_D.dYdY = diag(1./S2(:))-2*krontmp.*resprod;
            %ddJ_DdYdSigma
            J_D.dYdSigma = zeros([length(Y),length(Y)]);
            %ddJ_DdSigmadSigma
            J_D.dSigmadSigma = zeros([length(Y),length(Y)]);
            if nderiv >= 3
                %dddJ_DdYdYdY
                kron2tmp = bsxfun(@times,krontmp,permute(krontmp,[3,1,2]));
                tmp = -2*bsxfun(@times,(Y-Ym(:)).*N(:)./(S2(:).^2),permute(eye(length(Y)),[3,2,1])).*kron2tmp;
                tmp(isnan(tmp)) = 0;
                J_D.dYdYdY = 8*kron2tmp.*bsxfun(@times,resprod,permute(res,[3,2,1])) + tmp + permute(tmp,[3,1,2]) + permute(tmp,[2,3,1]);
                %dddJ_DdYdYdSigma
                J_D.dYdYdSigma = zeros([length(Y),length(Y),length(Y)]);
                %dddJ_DdYdSigmadSigma
                J_D.dYdSigmadSigma = zeros([length(Y),length(Y),length(Y)]);
                %dddJ_DdSigmadSigmadSigma
                J_D.dSigmadSigmadSigma = zeros([length(Y),length(Y),length(Y)]);
            end
        end
    end
end
end