function [SP,logL_m,logL_C,dlogL_mdxi,dlogL_Cdxi,ddlogL_mdxi2,ddlogL_Cdxi2] = logL_SCSH(xi, Model, Data, s, options)
% This routineshould work the same way as logL_PA does, except that from
% SCSH data, we are able to estimate the variance of the distribution
% (=variance coming from biological variability + measurement noise) by the
% SP approximation and compare it to the one computed from the SCSH data. 
% All we need is to estimate the variance of the biological variability
% along with the other parameters, or compute it from replicates.

    nderiv = 0.5 * (nargout-1) - 1;
    
    if isfield(Model.exp{s}, 'approx')
        optionsSim.approx = Model.exp{s}.approx;
    else
        optionsSim.approx = 'sp';
    end
    if isfield(Model.exp{s}, 'samples')
        optionsSim.samples = Model.exp{s}.samples;
    end
    
    %% Simulation
    if(nargout >= 4)
        [SP,my,Cy,dmydxi,dCydxi] = getSimulationSCSH(xi, Model, Data, s, optionsSim);
    else
        [SP,my,Cy] = getSimulationSCSH(xi, Model, Data, s, optionsSim);
    end
       
    %% Processing of simulation results, in the case that data points are missing, doubled, or multiple conditions are measured
    
    % Duplicate values in my if more than one data point at one time point
    if (size(Data{s}.condition,1) == 1)
        % If we do not have a dose response experiment
        if (size(Data{s}.SCSH.time,1) ~= size(my,1))
            k = 0;
            oldT = nan;
            tmp_my = nan(size(Data{s}.SCSH.time,1), size(my,2));
            tmp_Cy = nan(size(Data{s}.SCSH.time,1), size(Cy,2));
            if(nargout >= 4)
                tmp_dmydxi = nan(size(Data{s}.SCSH.time,1), size(my,2), size(dmydxi,3));
                tmp_dCydxi = nan(size(Data{s}.SCSH.time,1), size(dCydxi,2), size(dCydxi,3));
            end
            for j = 1 : size(Data{s}.SCSH.time,1)
                if (Data{s}.SCSH.time(j) ~= oldT)
                    k = k + 1;
                end
                tmp_my(j,:) = my(k,:);
                tmp_Cy(j,:) = Cy(k,:);
                if(nargout >= 4)
                    tmp_dmydxi(j,:,:) = dmydxi(k,:,:);
                    tmp_dCydxi(j,:,:) = dCydxi(k,:,:);
                end
                oldT = Data{s}.SCSH.time(j);
            end
            my = tmp_my;
            Cy = tmp_Cy;
            if(nargout >= 4)
                dmydxi = tmp_dmydxi;
                dCydxi = tmp_dCydxi;
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
    switch options.estimate_sigma
        case 0 
            % no estimation of noise parameters
            Sigma_m = Data{s}.SCSH.Sigma_m;
            Sigma_C = Data{s}.SCSH.C;
        case 1
            % standard estimation of noise parameters
            Sigma_m = Model.exp{s}.sigma_mean(Model.exp{s}.phi(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi)));
            Sigma_m = repmat(Sigma_m, [size(Data{s}.SCSH.m,1) 1]);
            Data{s}.SCSH.Sigma_m = Sigma_m;
            
            Sigma_C = Model.exp{s}.sigma_cov(Model.exp{s}.phi(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi)));
            Sigma_C = repmat(Sigma_C, [size(Data{s}.SCSH.C,1) 1]);
            Data{s}.SCSH.Sigma_C = Sigma_C;
            % Sigma_C = 0.1 * Data{s}.SCSH.C;
        case 2
            % optimal analytic computation of noise parameters
    end
  
    % Compute likelihood and derivatives for the mean
    switch Model.exp{s}.noise_model
        case 'normal'
            J_D_m = normal_noise(my(:), Data{s}.SCSH.m, Sigma_m, 1:size(Data{s}.SCSH.m, 1), nderiv);
        case 'lognormal'
            J_D_m = lognormal_noise(my(:), log(Data{s}.SCSH.m), Sigma_m, 1:size(Data{s}.SCSH.m, 1), nderiv);
    end

    % Compute likelihood and derivatives from biol. variability
    switch Model.exp{s}.variance_noise_model
        case 'normal'
            J_D_C = normal_noise(Cy(:), Data{s}.SCSH.C, Sigma_C, 1:size(Data{s}.SCSH.C, 1), nderiv);
        case 'lognormal'
            J_D_C = lognormal_noise(Cy(:), Data{s}.SCSH.C, Sigma_C, 1:size(Data{s}.SCSH.C, 1), nderiv);
    end

    % Write values to output
    logL_m = -J_D_m.val;
    if strcmp(optionsSim.approx, 'pa only')
        logL_C = 0;
    else
        logL_C = -J_D_C.val;
    end

    if (nderiv >= 1)
        % Compute derivative for dynamic parameters, scalings and offsets
        dlogL_mdy = reshape(-J_D_m.dY, size(Data{s}.SCSH.m));
        dlogL_mdxi = squeeze(nansum(nansum(repmat(dlogL_mdy, [1 1 size(dmydxi, 3)]) .* dmydxi, 2), 1));
        
        % Compute derivative for simga_mean parameters
        dlogL_mdSigma_m = reshape(-J_D_m.dSigma, size(Data{s}.SCSH.m));
        phi = Model.exp{s}.phi(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi));
        dSigma_mdphi = Model.exp{s}.dsigma_meandphi(phi);
        dSigma_mdbeta = permute(dSigma_mdphi, [2 3 1]) * Model.exp{s}.dphidbeta(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi));
        dSigma_mdxi = repmat(permute(dSigma_mdbeta * Model.exp{s}.dbetadxi(xi), [3 1 2]), [size(dlogL_mdSigma_m,1) 1 1]);
        dlogL_mdxi_Sigma_mPart = squeeze(nansum(nansum(repmat(dlogL_mdSigma_m, [1 1 length(xi)]) .* dSigma_mdxi, 2), 1));
        dlogL_mdxi = dlogL_mdxi + dlogL_mdxi_Sigma_mPart;
        
        if strcmp(optionsSim.approx, 'pa only')
            dlogL_Cdxi = zeros(size(xi));
        else
            % Compute derivative for dynamic parameters, scalings, offsets and
            % sigma_noise parameters
            dlogL_Cdy = reshape(-J_D_C.dY, size(Data{s}.SCSH.C));
            dlogL_Cdxi = squeeze(nansum(nansum(dCydxi .* repmat(dlogL_Cdy, [1 1 size(dmydxi, 3)]), 2), 1));

            % Compute derivative for simga_cov parameters
            dlogL_CdSigma_C = reshape(-J_D_C.dSigma, size(Data{s}.SCSH.C));
            dSigma_Cdphi = Model.exp{s}.dsigma_covdphi(phi);
            dSigma_Cdbeta = permute(dSigma_Cdphi, [2 3 1]) * Model.exp{s}.dphidbeta(Model.exp{s}.beta(xi), Model.exp{s}.delta(xi));
            dSigma_Cdxi = repmat(permute(dSigma_Cdbeta * Model.exp{s}.dbetadxi(xi), [3 1 2]), [size(dlogL_CdSigma_C,1) 1 1]);
            dlogL_Cdxi_Sigma_CPart = squeeze(nansum(nansum(repmat(dlogL_CdSigma_C, [1 1 length(xi)]) .* dSigma_Cdxi, 2), 1));
            dlogL_Cdxi = dlogL_Cdxi + dlogL_Cdxi_Sigma_CPart;
        end
        if (nderiv >= 2)
            switch Model.exp{s}.noise_model
                case 'normal'
                    % Term coming from sy' * sy
                    nan_ind = isnan(Data{s}.SCSH.m(:));
                    dres_mdxi = ((1 ./ Sigma_m(:)) * ones(1,length(xi))) .* reshape(dmydxi, numel(Data{s}.SCSH.m), length(xi));
                    dres_mdxi(nan_ind,:) = 0;
                    ddlogL_mdxi2 = -transpose(dres_mdxi) * dres_mdxi;
                    
                    if ~strcmp(optionsSim.approx, 'pa only')
                        % Term coming from sy' * sy
                        nan_ind = isnan(Data{s}.SCSH.C(:));
                        dres_Cdxi = ((1 ./ Sigma_C(:)) * ones(1,length(xi))) .* reshape(dCydxi, numel(Data{s}.SCSH.C), length(xi));
                        dres_Cdxi(nan_ind,:) = 0;
                        ddlogL_Cdxi2 = -transpose(dres_Cdxi) * dres_Cdxi;
                    end
                    
                    if (options.estimate_sigma == 1)
                        % Term 1 coming from s_sigma' * s_sigma
                        dresSigma_mdxi = ((1 ./ Sigma_m(:)) * ones(1,length(xi))) .* reshape(dSigma_mdxi, numel(Data{s}.SCSH.m), length(xi));
                        dresSigma_mdxi(nan_ind,:) = 0;
                        ddlogL_mdxi2 = ddlogL_mdxi2 + transpose(dresSigma_mdxi) * dresSigma_mdxi;

                        % Term 2 coming from s_sigma' * s_sigma
                        res = my(:) - Data{s}.SCSH.m(:);
                        dresSigma_mdxi2 = ((sqrt(3) * res ./ (Sigma_m(:).^2)) * ones(1,length(xi))) .* reshape(dSigma_mdxi, numel(Data{s}.SCSH.m), length(xi));
                        dresSigma_mdxi2(nan_ind,:) = 0;
                        ddlogL_mdxi2 = ddlogL_mdxi2 - transpose(dresSigma_mdxi2) * dresSigma_mdxi2;

                        if ~strcmp(optionsSim.approx, 'pa only')
                            % Term 1 coming from s_sigma' * s_sigma
                            dresSigma_Cdxi = ((1 ./ Sigma_C(:)) * ones(1,length(xi))) .* reshape(dSigma_Cdxi, numel(Data{s}.SCSH.C), length(xi));
                            dresSigma_Cdxi(nan_ind,:) = 0;
                            ddlogL_Cdxi2 = ddlogL_Cdxi2 + transpose(dresSigma_Cdxi) * dresSigma_Cdxi;

                            % Term 2 coming from s_sigma' * s_sigma
                            res = Cy(:) - Data{s}.SCSH.C(:);
                            dresSigma_Cdxi2 = ((sqrt(3) * res ./ (Sigma_C(:).^2)) * ones(1,length(xi))) .* reshape(dSigma_Cdxi, numel(Data{s}.SCSH.C), length(xi));
                            dresSigma_Cdxi2(nan_ind,:) = 0;
                            ddlogL_Cdxi2 = ddlogL_Cdxi2 - transpose(dresSigma_Cdxi2) * dresSigma_Cdxi2;
                        end
                    end

                case 'lognormal'
                    % To be done!
            end
        end
    end

    % Visulization
    if options.plot
        Sim_SCSH.m = my;
        Sim_SCSH.C = Cy;
        Sim_SCSH.Sigma_m = Sigma_m;
        Sim_SCSH.Sigma_C = Sigma_C;
        Sim_SCSH.t = Data{s}.SCSH.time;
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
