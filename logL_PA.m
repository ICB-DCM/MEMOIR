function [SP,logL_m,dlogL_mdxi,ddlogL_mdxi2] = logL_PA(xi, Model, Data, s, options)

    nderiv = nargout-2;
    
    if isfield(Model.exp{s}, 'approx')
        optionsSim.approx = Model.exp{s}.approx;
    else
        optionsSim.approx = 'sp';
    end
    if isfield(Model.exp{s}, 'samples')
        optionsSim.samples = Model.exp{s}.samples;
    end
    
    % Simulation
    if (nargout >= 3)
        [SP,my,dmydxi] = getSimulationPA(xi, Model, Data, s, optionsSim);
    else
        [SP,my] = getSimulationPA(xi, Model, Data, s, optionsSim);
    end

    %% Processing of simulation results, in the case that data points are missing, doubled, or multiple conditions are measured

    % Duplicate values in my if more than one data point at one time point
    if (size(Data{s}.condition,1) == 1)
        % No dose response experiment
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
        % Dose response experiment
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
    switch options.estimate_sigma
        case 0 
            % no estimation of noise parameters
            Sigma = Data{s}.PA.Sigma_m;
        case 1
            % standard estimation of noise parameters
            Sigma = Model.exp{s}.sigma_mean(Model.exp{s}.phi(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi)));
            Sigma = repmat(Sigma, [size(Data{s}.PA.m,1) 1]);
        case 2
            % optimal analytic computation of noise parameters
    end
    
    % Compute likelihood and derivatives for the mean
    switch Model.exp{s}.noise_model
        case 'normal'
            J_D = normal_noise(my(:), Data{s}.PA.m, Sigma, 1:size(Data{s}.PA.m, 1), min(nderiv, 1));
        case 'lognormal'
            J_D = lognormal_noise(my(:), Data{s}.PA.m, Sigma, 1:size(Data{s}.PA.m, 1), min(nderiv, 1));
    end
    
    % Write values to output
    logL_m = -J_D.val;
    if (nderiv >= 1)
        % Compute derivative for dynamic parameters, scalings and offsets
        dlogL_mdy = reshape(-J_D.dY, size(Data{s}.PA.m));
        dlogL_mdxi = squeeze(nansum(nansum(repmat(dlogL_mdy, [1 1 size(dmydxi, 3)]) .* dmydxi, 2), 1));
        
        % Compute derivative for simga_mean parameters
        dlogL_mdSigma = reshape(-J_D.dSigma, size(Data{s}.PA.m));
        phi = Model.exp{s}.phi(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi));
        dSigmadphi = Model.exp{s}.dsigma_meandphi(phi);
        dSigmadbeta = permute(dSigmadphi, [2 3 1]) * Model.exp{s}.dphidbeta(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi));
        dSigmadxi = repmat(permute(dSigmadbeta * Model.exp{s}.dbetadxi(xi), [3 1 2]), [size(dlogL_mdSigma,1) 1 1]);
        dlogL_mdxi_SigmaPart = squeeze(nansum(nansum(repmat(dlogL_mdSigma, [1 1 length(xi)]) .* dSigmadxi, 2), 1));
        dlogL_mdxi = dlogL_mdxi + dlogL_mdxi_SigmaPart;
        
        if (nderiv >= 2)
            switch Model.exp{s}.noise_model
                case 'normal'
                    % Term coming from sy' * sy
                    nan_ind = isnan(Data{s}.PA.m(:));
                    dres_mdxi = ((1 ./ Sigma(:)) * ones(1,length(xi))) .* reshape(dmydxi, numel(Data{s}.PA.m), length(xi));
                    dres_mdxi(nan_ind,:) = 0;
                    ddlogL_mdxi2 = -transpose(dres_mdxi) * dres_mdxi;
                    
                    % Term 1 coming from s_sigma' * s_sigma
                    dresSigmadxi = ((1 ./ Sigma(:)) * ones(1,length(xi))) .* reshape(dSigmadxi, numel(Data{s}.PA.m), length(xi));
                    dresSigmadxi(nan_ind,:) = 0;
                    ddlogL_mdxi2 = ddlogL_mdxi2 + transpose(dresSigmadxi) * dresSigmadxi;
                    
                    % Term 2 coming from s_sigma' * s_sigma
                    res = my(:) - Data{s}.PA.m(:);
                    dresSigmadxi2 = ((sqrt(3) * res ./ (Sigma(:).^2)) * ones(1,length(xi))) .* reshape(dSigmadxi, numel(Data{s}.PA.m), length(xi));
                    dresSigmadxi2(nan_ind,:) = 0;
                    ddlogL_mdxi2 = ddlogL_mdxi2 - transpose(dresSigmadxi2) * dresSigmadxi2;
                    
                case 'lognormal'
                    % To be done!
            end
        end
    end
    
    % Visualization
    if options.plot
        Sim_PA.m = my;
        Sim_PA.Sigma_m = Sigma;
        Sim_PA.t = Data{s}.PA.time;
        if isfield(Model.exp{s},'PA_post_processing_SP')
            Sim_PA.SP_max = SP.SP_max;
            Sim_PA.SP_min = SP.SP_min;
        else
            Sim_PA.SP_max = [];
            Sim_PA.SP_min = [];
        end
        Model.exp{s}.plot(Data{s}, Sim_PA, s);
    end
    
    % fprintf('\n\n %e', logL_m);
end
