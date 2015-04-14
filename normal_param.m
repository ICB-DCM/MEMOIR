%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR PARAMETER DENSITIES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_param(b,delta,type_D)

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,type_D);

if nargout >=1
    % J_b
    varargout{1} = 0.5*b'*invD*b ...
        + 0.5*log(det(D));
    if nargout >= 2
        % dJ_bdb
        varargout{2} = transpose(invD*b);
        % dJ_bddelta
        varargout{3} = transpose(0.5*squeeze(sum(sum(bsxfun(@times,dinvDddelta,bsxfun(@times,permute(b,[2,1]),permute(b,[1,2]))),1),2)) ... % 1/2*b*dinvD*b
            +0.5*squeeze(sum(sum(sum(bsxfun(@times,invD.*eye(length(b)),permute(dDddelta,[4,1,2,3])),2),1),3))); % 1/2*Tr(invD*dD)
        if nargout >= 4
            % ddJ_bdbdb
            varargout{4} = invD;
            % ddJ_bdbddelta
            varargout{5} = squeeze(sum(bsxfun(@times,dinvDddelta,b),2)); % dinvD*b it does not reall matter onto which of the first two dimensions we multiply b here as D and all its derivatives are symmetric
            % ddJ_bddeltaddelta
            varargout{6} = 0.5*squeeze(sum(sum(bsxfun(@times,ddinvDddeltaddelta,bsxfun(@times,permute(b,[2,1]),permute(b,[1,2]))),1),2)) ... % 1/2*b*ddinvD*b
                +0.5*squeeze(sum(sum(bsxfun(@times,squeeze(sum(bsxfun(@times,permute(dinvDddelta,[1,2,4,3]),permute(dDddelta,[4,1,2,5,3])),2)+sum(bsxfun(@times,invD,permute(ddDddeltaddelta,[5,1,2,3,4])),2)),eye(length(b))),1),2)); % 1/2*Tr(dinvD*dD + invD*ddD)
            if nargout >= 7
                % dddJ_bdbdbdb
                varargout{7} = zeros(length(b),length(b),length(b));
                % dddJ_bdbdbddelta
                varargout{8} = dinvDddelta;
                % dddJ_bdbddeltaddelta
                varargout{9} = squeeze(sum(bsxfun(@times,ddinvDddeltaddelta,b),2)); % dinvD*b it does not reall matter onto which of the first two dimensions we multiply b here as D and all its derivatives are symmetric
                if nargout >= 10
                    % ddddJ_bdbdbdbdb
                    varargout{10} = zeros(length(b),length(b),length(b),length(b));
                    % ddddJ_bdbdbdbddelta
                    varargout{11} = zeros(length(b),length(b),length(b),length(b));
                    % ddddJ_bdbdbddeltaddelta
                    varargout{12} = ddinvDddeltaddelta;
                    
                end
            end
        end
    end
end
end