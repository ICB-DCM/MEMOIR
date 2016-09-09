%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR PARAMETER DENSITIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_b = normal_param(b,delta,type_D,nderiv)

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);

if nderiv >= 0
    % J_b
    J_b.val = 0.5*b'*invD*b ...
        + 0.5*log(det(D));
    if nderiv >= 1
        % dJ_bdb
        J_b.db = transpose(invD*b);
        % dJ_bddelta
        J_b.ddelta = transpose(0.5*permute(sum(sum(bsxfun(@times,dinvDddelta,bsxfun(@times,permute(b,[2,1]),permute(b,[1,2]))),1),2),[1,3,2]) ... % 1/2*b*dinvD*b
            +0.5*permute(sum(sum(sum(bsxfun(@times,invD.*eye(length(b)),permute(dDddelta,[4,1,2,3])),2),1),3),[1,4,3,2])); % 1/2*Tr(invD*dD)
        if nderiv >= 2
            % ddJ_bdbdb
            J_b.dbdb = invD;
            % ddJ_bdbddelta
            J_b.dbddelta = permute(sum(bsxfun(@times,dinvDddelta,b),2),[1,3,2]); % dinvD*b it does not reall matter onto which of the first two dimensions we multiply b here as D and all its derivatives are symmetric
            % ddJ_bddeltaddelta
            J_b.ddeltaddelta = 0.5*permute(sum(sum(bsxfun(@times,ddinvDddeltaddelta,bsxfun(@times,permute(b,[2,1]),permute(b,[1,2]))),1),2),[3,4,1,2]) ... % 1/2*b*ddinvD*b
                +0.5*permute(sum(sum(bsxfun(@times,permute(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2),[1,3,4,5,2]),eye(length(b))),1),2),[3,4,1,2]); % 1/2*Tr(dinvD*dD + invD*ddD)
            if nderiv >= 3
                % dddJ_bdbdbdb
                J_b.dbdbdb = zeros(length(b),length(b),length(b));
                % dddJ_bdbdbddelta
                J_b.dbdbddelta = dinvDddelta;
                % dddJ_bdbddeltaddelta
                J_b.dbddeltaddelta = permute(sum(bsxfun(@times,ddinvDddeltaddelta,b),2),[1,3,4,2]); % dinvD*b it does not reall matter onto which of the first two dimensions we multiply b here as D and all its derivatives are symmetric
%                 if nargout >= 10
%                     % ddddJ_bdbdbdbdb
%                     varargout{10} = zeros(length(b),length(b),length(b),length(b));
%                     % ddddJ_bdbdbdbddelta
%                     varargout{11} = zeros(length(b),length(b),length(b),length(b));
%                     % ddddJ_bdbdbddeltaddelta
%                     varargout{12} = ddinvDddeltaddelta;
%                     
%                 end
            end
        end
    end
end
end