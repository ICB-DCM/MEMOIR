function [SP,logL_m,dlogL_mdxi,ddlogL_mdxi2] = logL_PA(xi, Model, Data, s, options)

    % Simulation
    if(nargout >= 3)
        [SP,my,dmydxi] = getSimulationPA(xi, Model, Data, s);
    else
        [SP,my] = getSimulationPA(xi, Model, Data, s);
    end

%     Log-Transformed Version
%     if(nargout >= 3)
%         for j = 1 : size(dmydxi,3)
%             dmydxi(:,:,j) = dmydxi(:,:,j) ./ (my(:,:) * log(10));
%         end
%     end
%     my = log10(my);



    %% Processing of simulation results, in the case that data points are missing, doubled, or multiple conditions are measured

    % Duplicate values in my if more than one data point at one time point
    if (size(Data{s}.condition,1) == 1)
        % If we do not have a dose response experiment
        if (size(Data{s}.PA.time,1) ~= size(my,1))
            k = 0;
            oldT = nan;
            tmp_my = nan(size(Data{s}.PA.time,1), size(my,2));
            if(nargout >= 3)
                tmp_dmydxi = nan(size(Data{s}.PA.time,1), size(my,2), size(dmydxi,3));
            end
            for j = 1 : size(Data{s}.PA.time,1)
                if (Data{s}.PA.time(j) ~= oldT), k = k + 1; end
                tmp_my(j,:) = my(k,:);
                if(nargout >= 3)
                    tmp_dmydxi(j,:,:) = dmydxi(k,:,:);
                end
                oldT = Data{s}.PA.time(j);
            end
            my = tmp_my;
            if(nargout >= 3)
                dmydxi = tmp_dmydxi;
            end
        end
    else
        % If we have a dose response experiment and multiple data points for one time point
        tmp_my = nan(size(Data{s}.condition,1), size(my,2));
        thisUniqueCondition = unique(Data{s}.condition, 'rows');
        if(nargout >= 3)
            tmp_dmydxi = nan(size(Data{s}.condition,1), size(my,2), size(dmydxi,3));
        end
        for j = 1 : size(Data{s}.condition,1) % number of conditions
            for iDose = 1 : size(thisUniqueCondition,1)
                if all(thisUniqueCondition(iDose,:)==Data{s}.condition(j,:))
                    tmp_my(j,:) =  my(iDose,:);
                    if(nargout >= 3)
                        tmp_dmydxi(j,:,:) = dmydxi(iDose,:,:);
                    end
                end
            end
        end
        my = tmp_my;
        if(nargout >= 3)
            dmydxi = tmp_dmydxi;
        end
    end
    
    
    
    %% Evaluation of the Likelihood
    
    % Evaluation of likelihood, likelihood gradient and hessian
    logL_m = - 0.5*nansum(nansum(((Data{s}.PA.m - my)./Data{s}.PA.Sigma_m).^2,1),2);
    fprintf('Nr: %2i,  LogL: %12.5f \n', s, logL_m);
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
        % Sim_PA.m = my;
        % Model.exp{s}.plot(Data{s},Sim_PA,s);
    end
    
end
