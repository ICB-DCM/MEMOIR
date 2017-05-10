function [SP,logL_m,logL_C,dlogL_mdxi,dlogL_Cdxi,ddlogL_mdxi2,ddlogL_Cdxi2] = logL_SCSH(xi, Model, Data, s, options)
% This routineshould work the same way as logL_PA does, except that from
% SCSH data, we are able to estimate the variance of the distribution
% (=variance coming from biological variability + measurement noise) by the
% SP approximation and compare it to the one computed from the SCSH data. 
% All we need is to estimate the variance of the biological variability
% along with the other parameters, or compute it from replicates.

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
    
    % Evaluation of likelihood, gradient and hessian from simulation
    logL_m = - 0.5*nansum(nansum(((Data{s}.SCSH.m - my)./Data{s}.SCSH.Sigma_m).^2,1),2);
    fprintf('Nr: %2i,  LogL: %12.5f \n', s, logL_m);
    if nargout >= 4
        dlogL_mdxi = squeeze(nansum(nansum(bsxfun(@times,(Data{s}.SCSH.m - my) ./ Data{s}.SCSH.Sigma_m.^2,dmydxi),1),2));
        
        if nargout >= 6
            wdmdxi = bsxfun(@times,1./Data{s}.SCSH.Sigma_m,SP.dmydxi);
            wdmdxi = reshape(wdmdxi,[numel(SP.my),size(SP.dmydxi,3)]);
            ddlogL_mdxi2 = -wdmdxi'*wdmdxi;
        end
    end
    
    % Evaluation of likelihood, gradient and hessian from biol. variability
    logL_C = - 0.5*nansum(nansum(nansum(((Data{s}.SCSH.C - SP.SCSH.Cy) ./ Data{s}.SCSH.Sigma_C).^2,1),2),3);
    fprintf('Nr: %2i,  LogL: %12.5f \n', s, logL_m);
    if nargout >= 4
        dlogL_Cdxi = squeeze(nansum(nansum(nansum(bsxfun(@times,(Data{s}.SCSH.C - Cy) ./ Data{s}.SCSH.Sigma_C.^2,dCydxi),1),2),3));
        
        if nargout >= 6
            % This will probably fail due to wrong dimensions...
            wdCdxi = bsxfun(@times,1./Data{s}.SCSH.Sigma_C,SP.dmydxi);
            wdCdxi = reshape(wdCdxi,[numel(SP.my),size(SP.dmydxi,3)]);
            ddlogL_Cdxi2 = -wdCdxi'*wdCdxi;
        end
    end
    
    % Visulization
    if options.plot
        % Sim_SCSH.m = my;
        % Sim_SCSH.C = Cy;
        % Model.exp{s}.plot(Data{s},Sim_SCSH,s);
    end
    
end
