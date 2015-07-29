%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR PARAMETER PENALTY   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = penal_param(b_s,delta,type_D,type_p,type_s)

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);

d_s = size(b_s,2);
n_b = size(b_s,1);

mu_s = 1/d_s*(sum(b_s,2));
S_s = (b_s*b_s');% + 1e-5*eye(n_b);

C_s = S_s - mu_s*mu_s';% + 1e-5*eye(n_b);
if type_p>3
    if type_s==1
        C_s = cov1para(C_s);
    elseif type_s==2
        C_s = shrinkDiag(C_s);
    end
end

VarC = 1/(d_s-1)*(diag(D)*diag(D)'+D.^2);

if nargout >=1
    
    if type_p==1
        % log(p({mu_s}|D))
        varargout{1} = - 0.5*log(det(D)) ... % log(det(D))
            - 0.5*d_s*(mu_s'*invD*mu_s); % mu*D^-1*mu
    elseif type_p==2||type_p==4
        % log(p({C_s}|D))
        varargout{1} = - 0.5*sum(log(VarC(:))) ...
            - 0.5*sum((C_s(:)-D(:)).^2./VarC(:));
    elseif type_p==3||type_p==5
        % log(p({mu_s}|D))
        part1 = - 0.5*log(det(D)) ... % log(det(D))
            - 0.5*d_s*(mu_s'*invD*mu_s); % mu*D^-1*mu
        % log(p({C_s}|D))
        part2 = - 0.5*sum(log(VarC(:))) ...
            - 0.5*sum((C_s(:)-D(:)).^2./VarC(:));
        varargout{1} = part1 + part2;
    end
    % J_s = log(p({mu_S,S_s}|D))
    %varargout{1} = + 0.5*(d_s-n_b-1)*log(det(S_s)) ... % log(det(S))
    %         - 0.5*(d_s+1)*log(det(D)) ... % log(det(D))
    %         - 0.5*d_s*(mu_s'*D*mu_s) ... % mu*D*mu
    %         - 0.5*trace(invD*S_s); % tr(invD*S)
    
    if nargout >= 2
        if type_p==1
            dmu_sdb_s = 1/d_s*repmat(eye(n_b),[1,1,d_s]);
            % dJ_sdb_s
            varargout{2} = - d_s*squeeze(sum(bsxfun(@times,mu_s'*invD,dmu_sdb_s),1));% - mu*invD*dmu;
            % dJ_sddelta
            varargout{3} = transpose(- 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dDddelta,[4,1,2,3])),2)),1),3)) ... % tr(invD*dD)
                - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,dinvDddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2))); % mu*dinvD*mu;
        elseif type_p==2||type_p==4
            dmu_sdb_s = 1/d_s*repmat(eye(n_b),[1,1,d_s]);
            dS_sdb_s = bsxfun(@times,permute(b_s,[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) + bsxfun(@times,permute(b_s,[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
            dCdb_s = dS_sdb_s - bsxfun(@times,permute(mu_s,[2,1]),permute(dmu_sdb_s,[1,4,2,3])) - bsxfun(@times,mu_s,permute(dmu_sdb_s,[4,1,2,3]));
            dVarCddelta = 1/(d_s-1)*bsxfun(@times,permute(diag(D),[3,2,1]),... % multiply with sigma^2 because delta = log(sigma^2)
                bsxfun(@times,bsxfun(@times,diag(D),permute(repmat(eye(n_b),[1,1,n_b]),[3 2 1])) + bsxfun(@times,diag(D)',permute(repmat(eye(n_b),[1,1,n_b]),[2 3 1])),...
                repmat(ones(n_b,n_b)+eye(n_b),[1,1,n_b]))); % multiply diagonal with 2
            % dJ_sdb_s
            varargout{2} =- squeeze(sum(sum(bsxfun(@times,repmat(bsxfun(@times,1./VarC,C_s - D),[1,1,n_b,d_s]),dCdb_s),1),2)); % - 1/VarC * (C_s - E(C_s)) * dCdb_s
            % dJ_sddelta
            varargout{3}  = transpose(squeeze(sum(sum(-.5*bsxfun(@times,1./VarC,dVarCddelta)... %-.5*1/VarC*dVarCddelta for case 'diag-matrix-logarithm'
                + bsxfun(@times,1./VarC,bsxfun(@times,C_s - D,dDddelta))...% + 1/VarC * (C_s - E(C_s) * dDddelta
                + bsxfun(@times,0.5*1./VarC.^2,bsxfun(@times,(C_s - D).^2,dVarCddelta)),1),2)));% + .5 * 1/VarC^2 * (C_s - E(C_s))^2 * dVarCddelta
            
        elseif type_p==3||type_p==5
            dmu_sdb_s = 1/d_s*repmat(eye(n_b),[1,1,d_s]);
            dS_sdb_s = bsxfun(@times,permute(b_s,[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) + bsxfun(@times,permute(b_s,[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
            dCdb_s = dS_sdb_s - bsxfun(@times,permute(mu_s,[2,1]),permute(dmu_sdb_s,[1,4,2,3])) - bsxfun(@times,mu_s,permute(dmu_sdb_s,[4,1,2,3]));
            dVarCddelta = 1/(d_s-1)*bsxfun(@times,permute(diag(D),[3,2,1]),... % multiply with sigma^2 because delta = log(sigma^2)
                bsxfun(@times,bsxfun(@times,diag(D),permute(repmat(eye(n_b),[1,1,n_b]),[3 2 1])) + bsxfun(@times,diag(D)',permute(repmat(eye(n_b),[1,1,n_b]),[2 3 1])),...
                repmat(ones(n_b,n_b)+eye(n_b),[1,1,n_b]))); % multiply diagonal with 2
            % dJ_sdb_s
            part1 = - d_s*squeeze(sum(bsxfun(@times,mu_s'*invD,dmu_sdb_s),1));% - mu*invD*dmu
            part2 = - squeeze(sum(sum(bsxfun(@times,repmat(bsxfun(@times,1./VarC,C_s - D),[1,1,n_b,d_s]),dCdb_s),1),2)); % - 1/VarC * (C_s - E(C_s)) * dCdb_s
            
            varargout{2} = part1 + part2;
            
            % dJ_sddelta
            part1 = transpose(- 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dDddelta,[4,1,2,3])),2)),1),3)) ... % tr(invD*dD)
                - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,dinvDddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2))); % mu*dinvD*mu
            
            part2 = transpose(squeeze(sum(sum(-.5*bsxfun(@times,1./VarC,dVarCddelta)... %-.5*1/VarC*dVarCddelta for case 'diag-matrix-logarithm'
                + bsxfun(@times,1./VarC,bsxfun(@times,C_s - D,dDddelta))...% + 1/VarC * (C_s - E(C_s) * dDddelta
                + bsxfun(@times,0.5*1./VarC.^2,bsxfun(@times,(C_s - D).^2,dVarCddelta)),1),2)));% + .5 * 1/VarC^2 * (C_s - E(C_s))^2 * dVarCddelta
            varargout{3} = part1 + part2;
        end
        
        %         invS_s = pinv(S_s);
        %         varargout{2} = 0.5*(d_s-n_b-1)*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invS_s,permute(dS_sdb_s,[5,1,2,3,4])),2)),1),3)) ... % tr(invS*dS)
        %             - d_s*squeeze(sum(bsxfun(@times,mu_s'*D,permute(dmu_sdb_s,[4,1,2,3])),2)) ... % mu*D*dmu
        %             - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dS_sdb_s,[5,1,2,3,4])),2)),1),3)); % tr(invD*dS)
        
        %         varargout{3} = transpose(-0.5*(d_s+1)*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(dDddelta,[4,1,2,3])),2)),1),3)) ... % tr(invD*dD)
        %             - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,dDddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)) ... % mu*dD*mu
        %             - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(S_s,[3,1,2]),permute(dinvDddelta,[1,2,4,3])),2)),1),3))); % tr(dinvD*S)
        
        if nargout >= 4
            
            % ddJ_sdb_sdb_s
            part1 = - d_s*squeeze(sum(sum(bsxfun(@times,bsxfun(@times,permute(dmu_sdb_s,[1,4,2,3]),invD),permute(dmu_sdb_s,[4,1,2,3])),2),1)); % dmu*invD*dmu
            % part2 = ;
            varargout{4} = part1;% + part2;
            
            % ddJ_sdb_sddelta
            part1 = - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,5,3]),bsxfun(@times,permute(mu_s,[2,1]),permute(dmu_sdb_s,[1,4,2,3]))),1),2)); % mu*dinvD*dmu
            % part2 = ;
            varargout{5} = part1;% + part2;
            
            % ddJ_sddeltaddelta
            part1 = -0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ...% 1/2*Tr(ddinvD*dD + invD*ddD)
                - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,ddinvDddeltaddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)); % mu*ddinvD*mu
            % part2 = ;
            varargout{6} = part1;% + part2;
            
            %             ddS_sdb_sdb_s = bsxfun(@times,permute(ones(size(b_s)),[1,3,4,2]),permute(eye(n_b),[3,1,2,4])) ...
            %                 + bsxfun(@times,permute(ones(size(b_s)),[3,1,4,2]),permute(eye(n_b),[1,3,2,4]));
            %             dinvS_sdb_s = -squeeze(sum(sum(bsxfun(@times,bsxfun(@times,invS_s,permute(dS_sdb_s,[5,1,2,6,3,4])),permute(invS_s,[3,4,1,2])),2),3));
            %             % ddJ_sdb_sdb_s
            %             varargout{4} = 0.5*(d_s-n_b-1)*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvS_sdb_s,[1,2,5,3,4]),permute(dS_sdb_s,[5,1,2,3,4])),2)+sum(bsxfun(@times,invS_s,permute(ddS_sdb_sdb_s,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ... % 1/2*Tr(dinvS*dS + invS*ddS)
            %                 - d_s*squeeze(sum(sum(bsxfun(@times,bsxfun(@times,permute(dmu_sdb_s,[1,4,2,3]),D),permute(dmu_sdb_s,[4,1,2,3])),2),1)) ...  % dmu*D*dmu
            %                 - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,invD,permute(ddS_sdb_sdb_s,[5,1,2,3,4])),2)),1),3)); % tr(invD*ddS)
            %             % ddJ_sdb_sddelta
            %             varargout{5} = - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,permute(dDddelta,[1,2,4,5,3]),bsxfun(@times,permute(mu_s,[2,1]),permute(dmu_sdb_s,[1,4,2,3]))),1),2)) ... % mu*dD*dmu
            %                 - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(dS_sdb_s,[1,2,5,3,4]),permute(dinvDddelta,[4,1,2,5,6,3])),2)),3),1)); % tr(dinvD*dS)
            %
            %             % ddJ_sddeltaddelta
            %             varargout{6} = -0.5*(d_s+1)*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2)),eye(n_b)),1),2)) ...% 1/2*Tr(ddinvD*dD + invD*ddD)
            %                 - 0.5*d_s*squeeze(sum(sum(bsxfun(@times,ddDddeltaddelta,bsxfun(@times,permute(mu_s,[2,1]),permute(mu_s,[1,2]))),1),2)) ... % mu*ddD*mu
            %                 - 0.5*squeeze(sum(sum(bsxfun(@times,permute(eye(n_b),[1,3,2]),sum(bsxfun(@times,permute(S_s,[3,1,2]),permute(ddinvDddeltaddelta,[1,2,5,3,4])),2)),1),3)); % tr(ddinvD*S)
        end
    end
end
end
