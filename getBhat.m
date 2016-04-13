function [B,G] = getBhat(xi, Model, Data, s, i, options, P_old)
    
    if(isfield(P_old{s},'SCTL'))
        if(isfield(P_old{s}.SCTL,'dbdxi'))
            bhat_si0 = P_old{s}.SCTL.bhat(:,i) + P_old{s}.SCTL.dbdxi(:,:,i)*(xi-xi_old);
        else
            bhat_si0 = P_old{s}.SCTL.bhat(:,i);
        end
    else
        bhat_si0 = zeros(length(Model.ind_b),1);
    end
    
    %% Estimation of single cell random effects
    % Higher order derivatives of the objective function for single cell parameters
    % here G is the Hessian of the objective function and bhat_si is the optimum of the objective function
    % with respect to b
    % F_diff and b_diff determine how many derivatives of the objective function and of the optimum need to
    % be computed
    
    % do multistart every few iterations
    fms = (mod(n_store,options.ms_iter)==0);
    
    switch(nderiv)
        case 0
            F_diff = 0;
            b_diff = 0;
            [bhat_si,G] ...
                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s,fms);
        case 1
            F_diff = 0;
            b_diff = 0;
            [bhat_si,G] ...
                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s,fms);
        case 2
            b_diff = 1;
            if(Model.integration)
                F_diff = 3;
                
                [bhat_si,dbhat_sidbeta,dbhat_siddelta,...
                    G,dGdb,pdGpdbeta,pdGpddelta] ...
                    = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s,fms);
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,2)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,3)
            else
                F_diff = 2;
                [bhat_si,dbhat_sidbeta,dbhat_siddelta,...
                    G] ...
                    = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s,fms);
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,2)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,3)
            end
        case 3
            if(Model.integration)
                F_diff = 4;
                b_diff = 2;
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,D,dDddelta,ddDddeltaddelta,invD,dinvDddelta,ddinvDddeltaddelta,t_s,Ym_si,ind,F_diff,b_diff,s),1e-4,2,4)
                [bhat_si,dbhat_sidbeta,dbhat_siddelta,ddbhat_sidbetadbeta,ddbhat_sidbetaddelta,ddbhat_siddeltaddelta,...
                    G,dGdb,pdGpdbeta,pdGpddelta,ddGdbdb,pddGdbpdbeta,pdpdGpdbetapdbeta,pddGdbpddelta,pdpdGpddeltapddelta,...
                    pdpdGpdbetapddelta] ...
                    = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s,fms);
            else
                F_diff = 3;
                b_diff = 2;
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,2)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,1,3)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(beta,@(beta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,2,4)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,2,5)
                % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(delta,@(delta)optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s),1e-4,3,6)
                [bhat_si,dbhat_sidbeta,dbhat_siddelta,ddbhat_sidbetadbeta,ddbhat_sidbetaddelta,ddbhat_siddeltaddelta,...
                    G,dGdb,pdGpdbeta] ...
                    = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s,fms);
            end
    end
    
    % Store bhat
    bhat(:,i) = bhat_si;
end