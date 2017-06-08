function [SP,logL_m,logL_C,dlogL_mdxi,dlogL_Cdxi,ddlogL_mdxi2,ddlogL_Cdxi2] = logL_SCSH(xi, Model, Data, s, options)
% This routineshould work the same way as logL_PA does, except that from
% SCSH data, we are able to estimate the variance of the distribution
% (=variance coming from biological variability + measurement noise) by the
% SP approximation and compare it to the one computed from the SCSH data. 
% All we need is to estimate the variance of the biological variability
% along with the other parameters, or compute it from replicates.

    nderiv = 0.5 * (nargout-1) - 1;

    %% Simulation
    if(nargout >= 4)
        [SP,my,Cy,dmydxi,dCydxi] = getSimulationSCSH(xi, Model, Data, s);
    else
        [SP,my,Cy] = getSimulationSCSH(xi, Model, Data, s);
    end
       
    %% Processing of simulation results, in the case that data points are missing, doubled, or multiple conditions are measured
    
    % Duplicate values in my if more than one data point at one time point
    if (size(Data{s}.condition,1) == 1)
        % If we do not have a dose response experiment
        if (size(Data{s}.SCSH.time,1) ~= size(my,1))
            k = 0;
            oldT = nan;
            tmp_my = nan(size(Data{s}.SCSH.time,1), size(my,2));
            if(nargout >= 4)
                tmp_dmydxi = nan(size(Data{s}.SCSH.time,1), size(my,2), size(dmydxi,3));
            end
            for j = 1 : size(Data{s}.SCSH.time,1)
                if (Data{s}.SCSH.time(j) ~= oldT)
                    k = k + 1;
                end
                tmp_my(j,:) = my(k,:);
                if(nargout >= 4)
                    tmp_dmydxi(j,:,:) = dmydxi(k,:,:);
                end
                oldT = Data{s}.SCSH.time(j);
            end
            my = tmp_my;
            if(nargout >= 4)
                dmydxi = tmp_dmydxi;
            end
        end
    else
        % If we have a dose response experiment and multiple data points for one time point
        tmp_my = nan(size(Data{s}.condition,1), size(my,2));
        thisUniqueCondition = unique(Data{s}.condition, 'rows');
        if(nargout >= 4)
            tmp_dmydxi = nan(size(Data{s}.condition,1), size(my,2), size(dmydxi,3));
        end
        for j = 1 : size(Data{s}.condition,1) % number of conditions
            for iDose = 1 : size(thisUniqueCondition,1)
                if all(thisUniqueCondition(iDose,:)==Data{s}.condition(j,:))
                    tmp_my(j,:) =  my(iDose,:);
                    if(nargout >= 4)
                        tmp_dmydxi(j,:,:) = dmydxi(iDose,:,:);
                    end
                end
            end
        end
        my = tmp_my;
        if(nargout >= 4)
            dmydxi = tmp_dmydxi;
        end
    end
    
    
    
    %% Evaluation of the Likelihood
    
    % Post-process data
    if isfield(Data{s}, 'SCSH_post_processing')
        [data_m, data_C] = Data{s}.SCSH_post_processing(Data{s}.SCSH.m, Data{s}.SCSH.C);
    end
    
    % Compute likelihood and derivatives for the mean
    switch Model.exp{s}.noise_model
        case 'normal'
            J_D_m = normal_noise(my(:), data_m, Data{s}.SCSH.Sigma_m, 1:size(data_m, 1), nderiv);
        case 'lognormal'
            J_D_m = lognormal_noise(my(:), log(data_m), Data{s}.SCSH.Sigma_m, 1:size(data_m, 1), nderiv);
    end
    
    % Compute likelihood and derivatives from biol. variability
    switch Model.exp{s}.variance_noise_model
        case 'normal'
            J_D_C = normal_noise(Cy(:), data_C, Data{s}.SCSH.Sigma_C, 1:size(data_m, 1), nderiv);
        case 'lognormal'
            J_D_C = lognormal_noise(Cy(:), data_C, Data{s}.SCSH.Sigma_C, 1:size(data_m, 1), nderiv);
    end
    
    % Write values to output
    logL_m = -J_D_m.val;
    logL_C = -J_D_C.val;
    if (logL_m == 0)
        logL_m = nan;
    end
    if (logL_C == 0)
        logL_C = nan;
    end
    if (nderiv >= 1)
        dlogL_mdy = reshape(-J_D_m.dY, size(data_m));
        dlogL_mdxi = squeeze(nansum(nansum(repmat(dlogL_mdy, [1 1 size(dmydxi, 3)]) .* dmydxi, 2), 1)); 
        dlogL_Cdy = reshape(-J_D_C.dY, size(data_C));
        dlogL_Cdxi = squeeze(nansum(nansum(nansum(dCydxi .* repmat(dlogL_Cdy, [1 1 1 size(dmydxi, 3)]), 3), 2), 1)); 
    end
    
%     % Evaluation of likelihood, gradient and hessian from simulation
%     logL_m = - 0.5*nansum(nansum(((data_m - my)./Data{s}.SCSH.Sigma_m).^2,1),2);
%     fprintf('Nr: %2i,  LogL: %12.5f \n', s, logL_m);
%     if nargout >= 4
%         dlogL_mdxi = squeeze(nansum(nansum(bsxfun(@times,(data_m - my) ./ Data{s}.SCSH.Sigma_m.^2,dmydxi),1),2));
%         
%         if nargout >= 6
%             wdmdxi = bsxfun(@times,1./Data{s}.SCSH.Sigma_m,SP.dmydxi);
%             wdmdxi = reshape(wdmdxi,[numel(SP.my),size(SP.dmydxi,3)]);
%             ddlogL_mdxi2 = -wdmdxi'*wdmdxi;
%         end
%     end
%     
%     % Evaluation of likelihood, gradient and hessian from biol. variability
%     logL_C = - 0.5*nansum(nansum(nansum(((Data{s}.SCSH.C - Cy) ./ Data{s}.SCSH.Sigma_C).^2,1),2),3);
%     fprintf('Nr: %2i,  LogL: %12.5f \n', s, logL_m);
%     if nargout >= 4
%         dlogL_Cdxi = squeeze(nansum(nansum(nansum(bsxfun(@times,(Data{s}.SCSH.C - Cy) ./ Data{s}.SCSH.Sigma_C.^2,dCydxi),1),2),3));
%         
%         if nargout >= 6
%             % This will probably fail due to wrong dimensions...
%             wdCdxi = bsxfun(@times,1./Data{s}.SCSH.Sigma_C,SP.dmydxi);
%             wdCdxi = reshape(wdCdxi,[numel(SP.my),size(SP.dmydxi,3)]);
%             ddlogL_Cdxi2 = -wdCdxi'*wdCdxi;
%         end
%     end

    % Visulization
    if options.plot
        Sim_SCSH.m = my;
        Sim_SCSH.C = Cy;
        if isfield(Model.exp{s},'SCSH_post_processing_SP')
            Sim_SCSH.SP_max = SP.SP_max;
            Sim_SCSH.SP_min = SP.SP_min;
        else
            Sim_SCSH.SP_max = [];
            Sim_SCSH.SP_min = [];
        end
        Model.exp{s}.plot(Data{s}, Sim_SCSH, s);
    end
    
end
