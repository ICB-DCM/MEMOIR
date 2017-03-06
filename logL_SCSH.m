function [P] = logL_SCSH(xi, Model, Data, s, options)
    
    error('currently not implemented (adapt from logLPA)')
    
    
    %% DONT FORGET TO ADD SIGMA_TECH TO CY
    
    % Covariance
    % for measurement noise we ignore the random effects at this point
    sigma = Model.exp{s}.sigma_noise(Model.exp{s}.phi(beta,zeros(n_b,1)));
    
    % adapt sigma to proper size
    if(op_SP.req(2))
        if(size(sigma,1) == n_t)
            if(size(sigma,2) == 1)
                C_tech = bsxfun(@times,repmat(sigma.^2,[1,n_y,n_y]),permute(eye(n_y),[3,1,2]));
            elseif(size(sigma,2) == n_y)
                C_tech = bsxfun(@times,repmat(sigma.^2,[1,1,n_y]),permute(eye(n_y),[3,1,2]));
            else
                error('Incompatible size of sigma parametrisation!')
            end
        elseif(size(sigma,2) == n_y)
            if(size(sigma,1) == 1)
                C_tech = bsxfun(@times,repmat(sigma.^2,[n_t,1,n_y]),permute(eye(n_y),[3,1,2]));
            else
                error('Incompatible size of sigma parametrisation!')
            end
        elseif(and(size(sigma,1)==1,size(sigma,2)==1))
            C_tech = bsxfun(@times,repmat(sigma.^2,[n_t,n_y,n_y]),permute(eye(n_y),[3,1,2]));
        else
            error('Incompatible size of sigma parametrisation!')
        end
    end
    
end

