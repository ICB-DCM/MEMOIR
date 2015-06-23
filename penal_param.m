%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR PARAMETER PENALTY   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function varargout = penal_param(b_s,delta,type_D)

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);

d_s = size(b_s,2);
n_b = size(b_s,1);

M = 1/d_s*sum(b_s,2);
dMdb_s = 1/d_s*repmat(eye(n_b),[1,1,d_s]);
C = cov(b_s',0);
dCdb_s = 1/(d_s-1)*bsxfun(@times,permute(b_s,[1,3,4,2]),permute(eye(n_b),[3,1,2,4])+permute(1/d_s*repmat(eye(n_b),[1,1,d_s]),[4,1,2,3])) ...
    + bsxfun(@times,permute(b_s,[3,1,4,2]),permute(eye(n_b),[1,3,2,4])+permute(1/d_s*repmat(eye(n_b),[1,1,d_s]),[1,4,2,3]));
dvarb_sdb_s = 1/(d_s-1)*2*(b_s-repmat(M,[1,d_s]));

Cx2y2 = var(b_s,1,2)*var(b_s,1,2)' + 2*cov(b_s',1);
Sigma_M = sqrt(var(b_s,0,2)/d_s) + 1e-5;

dSigma_Mdb_s
Sigma_C = sqrt(Cx2y2/ds - ((ds-2)*cov(b_s',1)-var(b_s,1,2)*var(b_s,1,2)')/(d_s*(d_s-1))) + 1e-5;

s_M = -log(Sigma_M.^2*2*pi);
s_C = -log(Sigma_C.^2*2*pi);
res_M =  - 0.5*((M)./Sigma_M).^2;
res_C =  - 0.5*((D-C)./Sigma_C).^2;

% d_s = size(b_s,2);
% n_b = size(b_s,1);
% 
% mu_s = 1/d_s*(sum(b_s,2));
% 
% S_s = b_s*(b_s)' + 1e-11*diag(ones(n_b,1));
% if nargout >=1
%     % J_s = -log(p({mu_S,S_s}|D))
%     varargout{1} = + 0.5*(d_s-n_b-1)*log(det(S_s)) ... % log(det(S))
%         - 0.5*(d_s+1)*log(det(D)) ... % log(det(D))
%         - 0.5*d_s*(mu_s'*D*mu_s) ... % mu*D*mu
%         - 0.5*trace(invD*S_s); % tr(invD*S)
%     
%     if nargout >= 2
%         dmu_sdb_s = 1/d_s*repmat(eye(n_b),[1,1,d_s]);
%         dS_sdb_s = bsxfun(@times,permute(b_s,[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) + bsxfun(@times,permute(b_s,[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
%         invS_s = pinv(S_s);
%         
%         % dJ_sdb_s
%         varargout{2} = 0.5*(d_s-n_b-1)*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invS_s,permute(dS_sdb_s,[5,1,2,3,4])),2)),1),3)) ... % tr(invS*dS)
%             - d_s*squeeze(sum(bsxfun(@times,mu_s'*D,permute(dmu_sdb_s,[4,1,2,3])),2)) ... % mu*D*dmu
%             - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dS_sdb_s,[5,1,2,3,4])),2)),1),3)); % tr(invD*dS)
%         
%         % dJ_sddelta
%         varargout{3} = transpose(-0.5*(d_s+1)*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dDddelta,[4,1,2,3])),2)),1),3)) ... % tr(invD*dD)
%             - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,dDddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)) ... % mu*dD*mu
%             - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(S_s,[3,1,2]),permute(dinvDddelta,[1,2,4,3])),2)),1),3))); % tr(dinvD*S)
%         if nargout >= 4
%             ddS_sdb_sdb_s = bsxfun(@times,permute(ones(size(b_s)),[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) ...
%                 + bsxfun(@times,permute(ones(size(b_s)),[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
%             dinvS_sdb_s = -squeeze(sum(sum(bsxfun(@times,bsxfun(@times,invS_s,permute(dS_sdb_s,[5,1,2,6,3,4])),permute(invS_s,[3,4,1,2])),2),3));
%             % ddJ_sdb_sdb_s
%             varargout{4} = 0.5*(d_s-n_b-1)*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvS_sdb_s,[1,2,5,3,4]),permute(dS_sdb_s,[5,1,2,3,4])),2)+sum(bsxfun(@times,invS_s,permute(ddS_sdb_sdb_s,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ... % 1/2*Tr(dinvS*dS + invS*ddS)
%                 - d_s*squeeze(sum(sum(bsxfun(@times,bsxfun(@times,permute(dmu_sdb_s,[1,4,2,3]),D),permute(dmu_sdb_s,[4,1,2,3])),2),1)) ...  % dmu*D*dmu
%                 - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(ddS_sdb_sdb_s,[5,1,2,3,4])),2)),1),3)); % tr(invD*ddS)
%             
%             
%             % ddJ_sdb_sddelta
%             varargout{5} = - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,permute(dDddelta,[1,2,4,5,3]),bsxfun(@times,permute(mu_s,[2,1]),permute(dmu_sdb_s,[1,4,2,3]))),1),2)) ... % mu*dD*dmu
%                 - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(dS_sdb_s,[1,2,5,3,4]),permute(dinvDddelta,[4,1,2,5,6,3])),2)),3),1)); % tr(dinvD*dS)
%             
%             % ddJ_sddeltaddelta
%             varargout{6} = -0.5*(d_s+1)*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ...% 1/2*Tr(ddinvD*dD + invD*ddD)
%                 - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,ddDddeltaddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)) ... % mu*ddD*mu
%                 - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(S_s,[3,1,2]),permute(ddinvDddeltaddelta,[1,2,5,3,4])),2)),1),3)); % tr(ddinvD*S)
%             
%         end
%     end
% end
end