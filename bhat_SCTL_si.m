% objective_SCTL_s1 is an auxiliary function for logL_CE_w_grad_2 and
% computes derivatives of the estimated random effect bhat wrt beta and
% delta the employed formulas can be derived via the implicit function
% theorem
%
% USAGE:
% ======
% [bhat,B.dbeta,B.ddelta,...] = bhat_SCTL_si(bhat,G,J.dbdbeta,J.dbddelta,~,~,~,G.db,G.pdbeta,G.pddelta,J.dbdbetadbeta,J.dbddeltaddelta,J.dbdbetaddelta)
%
% INPUTS:
% =======
% notation according to outputs of objective_SCTL_s1
%
% Outputs:
% ========
% estimated parameters bhat and derivatives wrt to beta and delta
%
% 2015/04/14 Fabian Froehlich

function B = bhat_SCTL_si(bhat,FIM,J,nderiv)
    
    invG = pinv(FIM.val);
    
    B.val = bhat;
    
    if nderiv>=1
        
        B.dbeta = -invG*squeeze(J.dbdbeta);
        B.ddelta = -invG*squeeze(J.dbddelta);
        
        if nderiv >= 2
            
            B.dbetadbeta = zeros(size(J.dbdbetadbeta));
            B.dbetaddelta = zeros(size(J.dbdbetaddelta));
            B.ddeltaddelta = zeros(size(J.dbddeltaddelta));
            
            for j=1:length(bhat)
                for k = 1:length(bhat)
                    if length(bhat) == 1
                        B.dbetadbeta = squeeze(B.dbetadbeta) - invG*(...
                            squeeze(J.dbdbetadbeta) ...
                            + transpose(B.dbeta)*transpose(squeeze(FIM.pdbeta(k,:,:))) ...
                            + transpose(transpose(B.dbeta)*transpose(squeeze(FIM.pdbeta(k,:,:)))) ...
                            + transpose(B.dbeta)*squeeze(FIM.db(k,:,:))*B.dbeta);
                        B.dbetadbeta = permute(B.dbetadbeta,[3,1,2]);
                        
                        B.dbetaddelta = transpose(squeeze(B.dbetaddelta)) - invG*(...
                            transpose(squeeze(J.dbdbetaddelta)) ...
                            + transpose(B.dbeta)*squeeze(FIM.pddelta(k,:,:)) ...
                            + transpose(B.ddelta)*squeeze(FIM.pdbeta(k,:,:)) ...
                            + transpose(B.dbeta)*squeeze(FIM.db(k,:,:))*B.ddelta);
                        B.dbetaddelta = permute(B.dbetaddelta,[3,1,2]);
                        
                        B.ddeltaddelta = squeeze(B.ddeltaddelta(j,:,:)) -invG(j,k)*(...
                            squeeze(J.dbddeltaddelta(k,:,:)) ...
                            + transpose(B.ddelta)*squeeze(FIM.pddelta(k,:,:)) ...
                            + transpose(transpose(B.ddelta)*squeeze(FIM.pddelta(k,:,:))) ...
                            + transpose(B.ddelta)*squeeze(FIM.db(k,:,:))*B.ddelta);
                        B.ddeltaddelta = permute(B.ddeltaddelta,[3,1,2]);
                    else
                        B.dbetadbeta(j,:,:) = squeeze(B.dbetadbeta(j,:,:)) - invG(j,k)*(...
                            squeeze(J.dbdbetadbeta(k,:,:)) ...
                            + transpose(B.dbeta)*squeeze(FIM.pdbeta(k,:,:)) ...
                            + transpose(transpose(B.dbeta)*squeeze(FIM.pdbeta(k,:,:))) ...
                            + transpose(B.dbeta)*squeeze(FIM.db(k,:,:))*B.dbeta);
                        
                        B.dbetaddelta(j,:,:) = squeeze(B.dbetaddelta(j,:,:)) - invG(j,k)*(...
                            squeeze(J.dbdbetaddelta(k,:,:)) ...
                            + transpose(B.dbeta)*squeeze(FIM.pddelta(k,:,:)) ...
                            + transpose(transpose(B.ddelta)*squeeze(FIM.pdbeta(k,:,:))) ...
                            + transpose(B.dbeta)*squeeze(FIM.db(k,:,:))*B.ddelta);
                        
                        B.ddeltaddelta(j,:,:) = squeeze(B.ddeltaddelta(j,:,:)) -invG(j,k)*(...
                            squeeze(J.dbddeltaddelta(k,:,:)) ...
                            + transpose(B.ddelta)*squeeze(FIM.pddelta(k,:,:)) ...
                            + transpose(transpose(B.ddelta)*squeeze(FIM.pddelta(k,:,:))) ...
                            + transpose(B.ddelta)*squeeze(FIM.db(k,:,:))*B.ddelta);
                    end
                end
            end
        end
    end
end