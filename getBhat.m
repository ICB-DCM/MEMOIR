function [B,G] = getBhat(xi, model, data, s, i, options, P_old)
    
    if(isfield(P_old{s},'SCTL'))
        if(isfield(P_old{s}.SCTL,'dbdxi'))
            bhat_si0 = P_old{s}.SCTL.bhat(:,i) + P_old{s}.SCTL.dbdxi(:,:,i)*(xi-P_old{s}.xi);
        else
            bhat_si0 = P_old{s}.SCTL.bhat(:,i);
        end
    else
        bhat_si0 = zeros(length(model.ind_b),1);
    end
    
    %% Estimation of single cell random effects
    % Higher order derivatives of the objective function for single cell parameters
    % here G is the Hessian of the objective function and bhat_si is the optimum of the objective function
    % with respect to b
    % F_diff and b_diff determine how many derivatives of the objective function and of the optimum need to
    % be computed
    
    beta = model.beta(xi);
    delta = model.delta(xi);
    
    b_diff = min(options.nderiv,2);
    if(options.nderiv>0)
        F_diff = 1+b_diff+options.integration;
    else
        F_diff = 0;
    end

    [B,G] = optimize_SCTL_si(model,data,bhat_si0,beta,delta,F_diff,b_diff,s,i,options,P_old);
    % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(model,data,bhat_si0,beta,D,dDddelta,ddDddeltaddelta,invD,dinvDddelta,ddinvDddeltaddelta,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,2,4)
end