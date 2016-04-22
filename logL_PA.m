function [SP,logL_m,dlogL_mdxi,ddlogL_mdxi2] = logL_PA(xi, Model, Data, s, options)

    % Simulation
    if(nargout >= 3)
        [SP,my,dmydxi] = getSimulationPA(xi, Model, Data, s);
    else
        [SP,my] = getSimulationPA(xi, Model, Data, s);
    end
    %         [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(theta) getSimulationPA(theta, Model, Data, s),1e-4,2,3);
    
    
    % Evaluation of likelihood, likelihood gradient and hessian
    logL_m = - 0.5*nansum(nansum(((Data{s}.PA.m - my)./Data{s}.PA.Sigma_m).^2,1),2);
    if nargout >= 3
        dlogL_mdxi = squeeze(nansum(nansum(bsxfun(@times,(Data{s}.PA.m - my)./Data{s}.PA.Sigma_m.^2,dmydxi),1),2));
        if nargout >= 4
            wdmdxi = bsxfun(@times,1./Data{s}.PA.Sigma_m,SP.dmydxi);
            wdmdxi = reshape(wdmdxi,[numel(SP.my),size(SP.dmydxi,3)]);
            ddlogL_mdxi2 = -wdmdxi'*wdmdxi;
        end
    end
    
    % Visulization
    if options.plot
        Sim_PA.m = my;
        Model.exp{s}.plot(Data{s},Sim_PA,s);
    end
    
end

