function [B,G] = getBhat(beta,delta, bhat_si0, model, data, s, i, options, P_old)
    
    %% Estimation of single cell random effects
    % Higher order derivatives of the objective function for single cell parameters
    % here G is the Hessian of the objective function and bhat_si is the optimum of the objective function
    % with respect to b
    % F_diff and b_diff determine how many derivatives of the objective function and of the optimum need to
    % be computed
    
    b_diff = min(options.nderiv,2);
    if(options.nderiv>0)
        F_diff = 1+b_diff+options.integration;
    else
        F_diff = 0;
    end

    [B,G] = optimize_SCTL_si(model,data,bhat_si0,beta,delta,F_diff,b_diff,s,i,options,P_old);
end