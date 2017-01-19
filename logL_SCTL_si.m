function [ logL, bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i)
%LOGL_SCTL_SI Summary of this function goes here
%   Detailed explanation goes here

if(isfield(P{s},'SCTL'))
    if(isfield(P{s}.SCTL,'dbdxi'))
        bhat_si0 = P{s}.SCTL.bhat(:,i) + P{s}.SCTL.dbdxi(:,:,i)*(xi-P{s}.xi);
    else
        bhat_si0 = P{s}.SCTL.bhat(:,i);
    end
else
    bhat_si0 = zeros(length(model.ind_b),1);
end

% compute bhat (optimum of the single cell likelihood
% (objective_SCTL_s1) with respect to b and the respective
% derivatives with respect to beta and delta
%
% testing:
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) getBhat_B( beta, delta, bhat_si0, model, data, s, i, options, P),1e-3,'val','dbeta')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) getBhat_B( beta, delta, bhat_si0, model, data, s, i, options, P),1e-3,'val','ddelta')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) getBhat_J( beta, delta, bhat_si0, model, data, s, i, options, P),1e-3,'val','dbeta')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) getBhat_J( beta, delta, bhat_si0, model, data, s, i, options, P),1e-3,'val','ddelta')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(beta,@(beta) getBhat_G( beta, delta, bhat_si0, model, data, s, i, options, P),1e-3,'val','dbeta')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(delta,@(delta) getBhat_G( beta, delta, bhat_si0, model, data, s, i, options, P),1e-3,'val','ddelta')

beta = model.beta(xi);
delta = model.delta(xi);
if(options.nderiv > 0)
    dbetadxi = model.dbetadxi(xi);
    ddeltadxi = model.ddeltadxi(xi);
    ddbetadxidxi = model.ddbetadxidxi(xi);
    dddeltadxidxi = model.dddeltadxidxi(xi);
end

[B,G,J,Sim] = getBhat( beta, delta, bhat_si0, model, data, s, i, options, P);

logL.D =  - J.J_D.val;
logL.T  = - J.J_T.val;
logL.b  = - J.J_b.val;
logL.val = -J.val;
bhat.val = B.val;

if(options.integration)
    % laplace approximation
    logL.I = - 0.5*log(abs(det(G.val)));
    logL.val = logL.val + logL.I;
end

if options.nderiv >= 1
    dbhat_sidbeta = B.dbeta;
    dbhat_siddelta = B.ddelta;
    % first order derivatives
    dbdxi = chainrule(dbhat_sidbeta,dbetadxi) + chainrule(dbhat_siddelta,ddeltadxi);
    bhat.dxi = dbdxi;

    logL.dxi = - chainrule(J.db,dbdxi) - chainrule(J.dbeta,dbetadxi) - chainrule(J.ddelta,ddeltadxi);
    
    if(options.integration)
        % laplace approximation
        invG = pinv(G.val);
        
        G.dbeta = G.dbeta + chainrule(G.db,dbhat_sidbeta);
        G.ddelta = G.ddelta + chainrule(G.db,dbhat_siddelta);
        G.dxi = chainrule(G.dbeta,dbetadxi) + chainrule(G.ddelta,ddeltadxi);
        logL.Idxi = - 0.5*permute(sum(sum(bsxfun(@times,permute(sum(bsxfun(@times,invG,permute(G.dxi,[4,1,2,3])),2),[1,3,4,2]),eye(length(bhat.val))),1),2),[1,3,2]); % 1/2*Tr(invG*dG)
        logL.dxi = logL.dxi + logL.Idxi; 
    end
%%  
%     if options.nderiv >= 2
%         % second order derivatives
%         
%         ddbhat_sidbetadbeta = B.dbetadbeta;
%         ddbhat_sidbetaddelta = B.dbetaddelta;
%         ddbhat_siddeltaddelta = B.ddeltaddelta;
%         
%         ddphidbdb = model.ddphidbdb(beta,bhat_si);
%         ddphidbdbeta = model.ddphidbdbeta(beta,bhat_si);
%         
%         ddphidbetadbeta = chainrule(dphidb,ddbhat_sidbetadbeta) + chainrule_ddxdydy_dydz(ddphidbdb,dbhat_sidbeta) ...
%             + permute(chainrule(permute(ddphidbdbeta,[1,3,2]),dbhat_sidbeta),[1,3,2]);
%         ddphidbetaddelta = chainrule(dphidb,ddbhat_sidbetaddelta) + chainrule_ddxdydy_dydz_dydv(ddphidbdb,dbhat_sidbeta,dbhat_siddelta);
%         ddphiddeltaddelta = chainrule(dphidb,ddbhat_siddeltaddelta) + chainrule_ddxdydy_dydz(ddphidbdb,dbhat_siddelta);
%         
%         if(numel(bhat_si)==1) % we have to do this manually here since 3rd order tensors with trailing 1 dimensional orders are not possible in matlab ...
%             ddphidxidxi = chainrule(dphidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddphidbetadbeta,dbetadxi) ...
%                 + chainrule(dphiddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddphiddeltaddelta,ddeltadxi) ...
%                 + permute(sum(bsxfun(@times,bsxfun(@times,ddphidbetaddelta,permute(dbetadxi,[3,1,2])),permute(ddeltadxi,[3,4,5,2,1])),2),[1,3,4,2]) ...
%                 + permute(sum(bsxfun(@times,bsxfun(@times,ddphidbetaddelta,permute(dbetadxi,[3,1,2])),permute(ddeltadxi,[3,4,5,2,1])),2),[1,4,3,2]);
%         else
%             ddphidxidxi = chainrule(dphidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddphidbetadbeta,dbetadxi) ...
%                 + chainrule(dphiddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddphiddeltaddelta,ddeltadxi) ...
%                 + chainrule_ddxdydy_dydz_dydv(ddphidbetaddelta,dbetadxi,ddeltadxi) + permute(chainrule_ddxdydy_dydz_dydv(ddphidbetaddelta,dbetadxi,ddeltadxi),[1,3,2]);
%         end
%         
%         ddJ_DdphidY = bsxfun(@times,ddJ_DdYdY,permute(dY_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_DdYdSigma,permute(dSigma_noisedphi,[3,1,2]));
%         
%         ddJ_DdphidSigma = bsxfun(@times,ddJ_DdYdSigma,permute(dY_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_DdSigmadSigma,permute(dSigma_noisedphi,[3,1,2]));
%         
%         ddJ_Ddphidphi = chainrule(dJ_DdY,ddY_sidphidphi) ...
%             + squeeze(sum(bsxfun(@times,ddJ_DdphidY,permute(dY_sidphi,[3,1,4,2])) ...
%             + bsxfun(@times,ddJ_DdphidSigma,permute(dSigma_noisedphi,[3,1,4,2])),2));
%         
%         ddJ_Ddxidxi = chainrule(dJ_Ddphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Ddphidphi,dphidxi);
%         
%         ddJ_TdphidT = bsxfun(@times,ddJ_TdTdT,permute(dT_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_TdTdR,permute(dR_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_TdTdSigma,permute(dSigma_timedphi,[3,1,2]));
%         
%         ddJ_TdphidR = bsxfun(@times,ddJ_TdTdR,permute(dR_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_TdRdR,permute(dR_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_TdRdSigma,permute(dSigma_timedphi,[3,1,2]));
%         
%         ddJ_TdphidSigma = bsxfun(@times,ddJ_TdTdSigma,permute(dT_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_TdRdSigma,permute(dR_sidphi,[3,1,2])) + ...
%             bsxfun(@times,ddJ_TdSigmadSigma,permute(dSigma_timedphi,[3,1,2]));
%         
%         ddJ_Tdphidphi = chainrule(dJ_TdT,ddT_sidphidphi) ...
%             + chainrule(dJ_TdR,ddR_sidphidphi) ...
%             + squeeze(sum(bsxfun(@times,ddJ_TdphidT,permute(dT_sidphi,[3,1,4,2])) ...
%             + bsxfun(@times,ddJ_TdphidR,permute(dR_sidphi,[3,1,4,2])) ...
%             + bsxfun(@times,ddJ_TdphidSigma,permute(dSigma_timedphi,[3,1,4,2])),2));
%         
%         ddJ_Tdxidxi = chainrule(dJ_Tdphi,ddphidxidxi) + chainrule_ddxdydy_dydz(ddJ_Tdphidphi,dphidxi);
%         
%         ddbdxidxi = chainrule(dbhat_sidbeta,ddbetadxidxi) + chainrule_ddxdydy_dydz(ddbhat_sidbetadbeta,dbetadxi) ...
%             + chainrule(dbhat_siddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(ddbhat_siddeltaddelta,ddeltadxi) ...
%             + chainrule_ddxdydy_dydz_dydv(ddbhat_sidbetaddelta,dbetadxi,ddeltadxi) ...
%             + chainrule_ddxdydy_dydz_dydv(permute(ddbhat_sidbetaddelta,[1,3,2]),ddeltadxi,dbetadxi);
%         
%         bhat.dxidxi = ddbdxidxi;
%         
%         ddJ_bdxidxi = chainrule(dJ_bdb,ddbdxidxi) + chainrule_ddxdydy_dydz(ddJ_bdbdb,dbdxi) ...
%             + chainrule(pdJ_bpddelta,dddeltadxidxi) + chainrule_ddxdydy_dydz(pdpdJ_bpddeltapddelta,ddeltadxi) ...
%             + chainrule_ddxdydy_dydz_dydv(dpdJ_bdbpddelta,dbdxi,ddeltadxi) ...
%             + chainrule_ddxdydy_dydz_dydv(permute(dpdJ_bdbpddelta,[2,1]),ddeltadxi,dbdxi);
%         
%         logL_D.dxidxi = - ddJ_Ddxidxi;
%         logL_T.dxidxi = - ddJ_Tdxidxi;
%         logL_b.dxidxi = - ddJ_bdxidxi;
%         
%         if(options.integration)
%             % laplace approximation
%             invG = pinv(G.val);
%             
%             ddGdbetadbeta = G.dbetadbeta + 2*chainrule(permute(G.dbdbeta,[1,2,4,3]),dbhat_sidbeta) ...
%                 + chainrule_ddxdydy_dydz(G.dbdb,dbhat_sidbeta) + chainrule(G.db,ddbhat_sidbetadbeta);
%             
%             ddGddeltaddelta = G.ddeltaddelta + 2*chainrule(permute(G.dbddelta,[1,2,4,3]),dbhat_siddelta) ...
%                 + chainrule_ddxdydy_dydz(G.dbdb,dbhat_siddelta) + chainrule(G.db,ddbhat_siddeltaddelta);
%             
%             ddGdbetaddelta = G.dbetaddelta + chainrule(permute(G.dbdbeta,[1,2,4,3]),dbhat_siddelta) ...
%                 + permute(chainrule(permute(G.dbddelta,[1,2,4,3]),dbhat_sidbeta),[1,2,4,3]) ...
%                 + chainrule_ddxdydy_dydz_dydv(G.dbdb,dbhat_sidbeta,dbhat_siddelta) ...
%                 + chainrule(G.db,ddbhat_sidbetaddelta);
%             
%             ddGdxidxi = chainrule_ddxdydy_dydz(ddGdbetadbeta,dbetadxi) + chainrule(G.dbeta,ddbetadxidxi) ...
%                 + chainrule_ddxdydy_dydz(ddGddeltaddelta,ddeltadxi) + chainrule(dGddelta,dddeltadxidxi) ...
%                 + 2*chainrule_ddxdydy_dydz_dydv(ddGdbetaddelta,dbetadxi,ddeltadxi);
%             
%             dinvGdxi = squeeze(sum(bsxfun(@times,invG,permute(squeeze(sum(bsxfun(@times,permute(dGdxi,[4,1,2,3]),invG),2)),[4,1,2,3])),2));
%             
%             logL_I.dxidxi = - 0.5*squeeze(sum(sum(squeeze(bsxfun(@times,sum(bsxfun(@times,permute(dinvGdxi,[1,2,4,3]),permute(dGdxi,[4,1,2,5,3])),2)+sum(bsxfun(@times,invG,permute(ddGdxidxi,[5,1,2,3,4])),2),permute(eye(length(bhat_si)),[1,3,2]))),1),2)); % 1/2*Tr(dinvG*dg + invG*ddG)
%         end
%         
%     end
%%    
end
end

function B = getBhat_B( beta, delta, bhat_si0, model, data, s, i, options, P)
    [B,G,J,Sim] = getBhat( beta, delta, bhat_si0, model, data, s, i, options, P);
end

function G = getBhat_G( beta, delta, bhat_si0, model, data, s, i, options, P)
    [B,G,J,Sim] = getBhat( beta, delta, bhat_si0, model, data, s, i, options, P);
end

function J = getBhat_J( beta, delta, bhat_si0, model, data, s, i, options, P)
    [B,G,J,Sim] = getBhat( beta, delta, bhat_si0, model, data, s, i, options, P);
end


