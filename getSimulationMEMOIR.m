function getSimulationMEMOIR(xi, Model, Data, ind_exp, ind_xi, options)
% This routine plots a number of specified experiments (ind_exp) with a 
% fine time grid and analyzes relative sensitivities for a given parameter 
% set (ind_xi)

    %% Preprocess inputs
    if (~isfield(options, 'sensi'))
        options.sensi = 0;
    end
    if (~isfield(options, 'fine'))
        options.fine = 1;
    end
    if (~isfield(options, 'approx'))
        options.approx = 'sp';
    end
    if (~isfield(options, 'nSamples'))
        options.nSamples = 1000;
    end
    if (~isfield(options, 'samples'))
        options.samples = [];
    end
    
%     if (options.sensi > 0)
%         % Collect all measurands
%         measurands = cell(1,0);
%         for s = ind_exp
%             measurands = {measurands{1:size(measurands,2)}, Data{s}.measurands{1:length(Data{s}.measurands)}};
%         end
%         measurands = unique(measurands);
%         
%         % Get length of parameter vector
%         parameters = 1 : length(Model.param);
%     
%         % Check validity of parameters for sensitivities to be computed
%         for j = 1 : length(ind_xi)
%             if ~any(parameters == ind_xi(j))
%                 error('The vector of parameter indices for sensitivity analysis is invalid.');
%             end
%         end
%         
%         % Check validity of observables for sensitivities to be computed
%         for iMeas = options.measurands
%             if ~any(strcmp(measurands, iMeas))
%                 err_msg = ['The observable '  iMeas{1,1} ...
%                     ' for which sensitivity is requested can not be found among the observables of the specified experiments.'];
%                 error(err_msg);
%             end
%         end
%         
%         %intSensiExp
%     end
    

    %% Perform run through all experiments
    % long_Y_true = nan(8,10,6,10000);
    
    for s = ind_exp
    % --- Loop over experiments -------------------------------------------
 
        % Get data type
        if isfield(Data{s}, 'SCTLstat')
            exp_type = 'SCTLstat';
        elseif isfield(Data{s}, 'SCSH')
            exp_type = 'SCSH';
            req = [1,1,0,0,0,1,0];
        elseif isfield(Data{s}, 'PA')
            exp_type = 'PA';
            req = [1,0,0,0,0,1,0];
        end
        post_proc = [exp_type '_post_processing'];
        
        % Set time vector accordingly
        if options.fine
            t_sim = linspace(Data{s}.(exp_type).time(1), Data{s}.(exp_type).time(end), 1000)';
            t_sim = unique([t_sim; Data{s}.(exp_type).time]);
            t_ind = [];
            for iT = 1 : length(t_sim)
                if any(t_sim(iT) == Data{s}.(exp_type).time)
                    t_ind = [t_ind, iT];
                end
            end
        else
            t_sim = Data{s}.(exp_type).time;
            t_ind = 1 : length(unique(Data{s}.(exp_type).time));
        end

        conditions = unique(Data{s}.condition,'rows');
        if (size(conditions,1) > 1)
            t_ind = 1:size(conditions,1);
        end
        
        % Preallocate arrays with simulation results
        my = [];
        Cy = [];
        my_true = [];
        Cy_true = [];
        
        for iCondition = 1 : size(conditions,1)
            %% Simulation of the model
            % Set options for sigma point routine
            op_SP.nderiv = options.sensi;
            op_SP.req = req;
            op_SP.type_D = Model.type_D;
            op_SP.approx = 'sp';
            op_SP.plot = 0;
            op_SP.nsamples = options.nSamples;

            if strcmp(options.approx, 'sp')
                % Call simulation
                SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model, t_sim, phi, conditions(iCondition,:), Model.exp{s}.scale), ...  = nonfun (in getSigmaPointApp)
                    xi, ...
                    Model.exp{s}, ... = estruct (in getSigmaPointApp)
                    op_SP);
            elseif strcmp(options.approx, 'samples')
                % Call simulation
                op_SP.approx = 'samples';
                if ~isempty(options.samples)
                    op_SP.samples = options.samples;
                end
                SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model, t_sim, phi, conditions(iCondition,:), Model.exp{s}.scale), ...  = nonfun (in getSigmaPointApp)
                    xi, ...
                    Model.exp{s}, ... = estruct (in getSigmaPointApp)
                    op_SP);
                 SP.my_true = SP.my;
                 SP.Y_true = SP.Y;
                 if isfield(SP, 'Cy')
                    SP.Cy_true = SP.Cy;
                 end
            elseif strcmp(options.approx, 'both')
                % Call simulation
                SP = testSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model, t_sim, phi, conditions(iCondition,:), Model.exp{s}.scale), ...  = nonfun (in getSigmaPointApp)
                    xi, ...
                    Model.exp{s}, ... = estruct (in getSigmaPointApp)
                    op_SP);
            else
                % Call simulation
                op_SP.approx = Model.exp{s}.approx;
                if ~isempty(Model.exp{s}.samples)
                    op_SP.samples = Model.exp{s}.samples;
                end
                SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model, t_sim, phi, conditions(iCondition,:), Model.exp{s}.scale), ...  = nonfun (in getSigmaPointApp)
                    xi, ...
                    Model.exp{s}, ... = estruct (in getSigmaPointApp)
                    op_SP);
            end

            %% Post-processing
            % Store the simulation results and apply scaling
            switch Model.exp{s}.scale
                case 'log'
                    % TBD!
                    tmp = arrayfun(@(x) diag(squeeze(SP.my(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                    my = [my; exp(SP.my + transpose([tmp{:}])/2)];
                    if ~strcmp(options.approx, 'sp')
                        tmp_true = arrayfun(@(x) diag(squeeze(SP.my_true(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                        my_true = [my_true; exp(SP.my_true + transpose([tmp_true{:}])/2)];
                    end

                    if strcmp(exp_type, 'SCSH')
                        tmpCy = bsxfun(@plus, repmat(tmp, 1, 1, size(SP.Cy,3)), permute(repmat(tmp, 1, 1, size(SP.Cy,3)), [1,3,2]));
                        Cy = [Cy; exp(tmpCy) .* (exp(SP.Cy) - ones(size(SP.Cy)))];
                        if ~strcmp(options.approx, 'sp')
                            tmpCy_true = bsxfun(@plus, repmat(tmp_true, 1, 1, size(SP.Cy_true,3)), permute(repmat(tmp_true, 1, 1, size(SP.Cy_true,3)), [1,3,2]));
                            Cy_true = [Cy_true; exp(tmpCy_true) .*  (exp(SP.Cy_true) - ones(size(SP.Cy_true)))];
                        end
                    end

                case 'log10'
                    % TBD!
                    tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                    my = [my; 10.^(SP.my + transpose([tmp{:}])/2)];
                    if strcmp(exp_type, 'SCSH')
                        tmpCy = bsxfun(@plus, repmat(tmp, 1, 1, size(SP.Cy,3)), permute(repmat(tmp, 1, 1, size(SP.Cy,3)), [1,3,2]));
                        Cy = [Cy; 10.^(tmpCy) .*  (10.^(SP.Cy) - ones(size(SP.Cy)))];
                    end

                case 'lin'
                    my = [my; SP.my];
                    if any(strcmp(options.approx, {'samples', 'both'}))
                        my_true = [my_true; SP.my_true];
                    end
                    if strcmp(exp_type, 'SCSH')
                        Cy = [Cy; SP.Cy];
                        if any(strcmp(options.approx, {'samples', 'both'}))
                            Cy_true = [Cy_true; SP.Cy_true];
                        end
                    end
            end
        end

        % Store gradients of means and variances
        if(options.sensi > 0)
            switch Model.exp{s}.scale
                case 'log'
                    % Set sizes for arrays
                    nt = size(SP.dmydxi,1);
                    np = size(SP.dmydxi,4);
                    ny = size(SP.dmydxi,2);

                    % To be checked!
                    dtmpdxi = arrayfun(@(x,y) diag(squeeze(SP.dCydxi(x,:,:,y))),repmat(1:nt,[np,1]),...
                        repmat(transpose(1:np),[1,nt]),'UniformOutput',false);
                    dmydxi = bsxfun(@times,my,SP.dmydxi) ...
                        + bsxfun(@times,my,permute(reshape([dtmpdxi{:}]/2,...
                        [ny,np,nt]),[3,1,2]));
                    dmydxi(isnan(dmydxi)) = 0;
                    if strcmp(exp_type, 'SCSH')
                    end

                case 'log10'
                    % To be checked!
                    tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                    my = 10.^(SP.my + transpose([tmp{:}])/2);
                    if strcmp(exp_type, 'SCSH')
                    end

                case 'lin'
                    dmydxi = SP.dmydxi;
                    if strcmp(exp_type, 'SCSH')
                        dCydxi = SP.dCydxi;
                    end
            end
        else
            dmydxi = zeros([size(my,1),size(my,2),length(xi)]);
            if strcmp(exp_type, 'SCSH')
                dCydxi = zeros([size(my,1),size(my,2),size(my,2),length(xi)]);
            end
        end

        % Apply additional user-defined post-processing
        if isfield(Model.exp{s}, post_proc)
            if(options.sensi > 0)
                SP.dmydxi = zeros([size(SP.my) size(xi,1)]);
            end
            if strcmp(exp_type, 'PA')
                [my,dmydxi] = feval(Model.exp{s}.(post_proc), my, dmydxi, xi);
                if any(strcmp(options.approx, {'samples', 'both'}))
                    %[my_true,~] = feval(Model.exp{s}.(post_proc), my_true, dmydxi);
                end
            elseif strcmp(exp_type, 'SCSH')
                [my,Cy,dmydxi,dCydxi] = feval(Model.exp{s}.(post_proc), my, Cy, dmydxi, dCydxi);
                if any(strcmp(options.approx, {'samples', 'both'}))
                    [my_true,Cy_true,~,~] = feval(Model.exp{s}.(post_proc), my_true, Cy_true, dmydxi, dCydxi);
                end
            end
        end


        %% Do sensitivity analysis

        if (options.sensi > 0)
            % Find indices of measurands for this experiment
            measInd = [];
            for iMeas = options.measurands
                for j = 1 : length(Data{s}.measurands)
                    if strcmp(Data{s}.measurands{j}, iMeas)
                        measInd = [measInd, j];
                    end
                end
            end

            % Compute relative sensitivities by (very simple) intgeration
            for iT = 1:size(dCydxi,1)
                tmp = permute(dCydxi(iT,1,1,:), [4 1 2 3]);
                tmp_dCydxi(iT,1,:) = permute(tmp, [2 3 1]);
            end
            dAlldxi = [dmydxi, tmp_dCydxi];
            relSensiExp = dAlldxi(:,measInd,ind_xi);% ./ repmat(my, [1 1 length(ind_xi)]);
            intSensiExp = zeros(size(relSensiExp,2),size(relSensiExp,3));
            for j = 1 : size(relSensiExp,1)
                intSensiExp = intSensiExp + (Sim.t(j+1)-Sim.t(j)) * 0.5 * permute((relSensiExp(j,:,:)+relSensiExp(j+1,:,:)), [2 3 1]);
            end
        end

        %% Plotting of fine simulation
        % Assign values for plotting
        if any(strcmp(options.approx, {'samples', 'both'}))
            Sim.mFineTrue = my_true;
        end
        Sim.mFine = my;
        Sim.t = unique(t_sim);
        SigmaStruct = processSigma(Data{s}, Sim.mFine, [], exp_type);
        Sim.Sigma_m = SigmaStruct.Sigma_m;
        if strcmp(exp_type, 'SCSH')
            if any(strcmp(options.approx, {'samples', 'both'}))
                Sim.CFineTrue = Cy_true;
                Sim.Y_true = SP.Y_true;
            end
            Sim.CFine = Cy;
            SigmaStruct = processSigma(Data{s}, Sim.mFine, Cy, exp_type);
            Sim.Sigma_C = SigmaStruct.Sigma_C;
        end

        % Process values in Sim.m for residual plot
        Sim = processSimulation(Sim, t_ind, Data{s}, exp_type);

        % long_Y_true(s-21,:,:,:) = Sim.Y_true;
        % Plotting
        Model.exp{s}.plot(Data{s}, Sim, s);  

        % Clean up
        clear Sim;
    % --- End of loop over experiments ------------------------------------
    end
    % save('simData.mat', 'long_Y_true');
end



function SigmaStruct = processSigma(thisData, mdata, Cdata, type)

    % Mean sigma
    sigma_m = thisData.(type).Sigma_m(1,:);
    
    if(size(sigma_m,1) == size(mdata,1))
        if(size(sigma_m,2) == 1)
            Sigma_m = repmat(sigma_m,[1,size(mdata,2)]);
        elseif(size(sigma_m,2) == size(mdata,2))
            Sigma_m = sigma_m;
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma_m,2) == size(mdata,2))
        if(size(sigma_m,1) == 1)
            Sigma_m = repmat(sigma_m,[size(mdata,1),1]);
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma_m,1)==1,size(sigma_m,2)==1))
        Sigma_m = repmat(sigma_m,size(mdata));
    else
        error('Incompatible size of sigma parametrisation!')
    end
    SigmaStruct.Sigma_m = Sigma_m;
        
    if (strcmp(type, 'SCSH') && ~isempty(Cdata))
        % Single-Cell snapshotsigma
        sigma_C = thisData.(type).Sigma_C(1,:,:);
        
        if(size(sigma_C,1) == size(Cdata,1))
            if(size(sigma_C,2) == 1)
                Sigma_C = repmat(sigma_C,[1,size(Cdata,2)]);
            elseif(size(sigma_C,2) == size(Cdata,2))
                Sigma_C = sigma_C;
            else
                error('Incompatible size of sigma parametrisation!')
            end
        elseif(size(sigma_C,2) == size(Cdata,2))
            if(size(sigma_C,1) == 1)
                Sigma_C = repmat(sigma_C,[size(Cdata,1),1]);
            else
                error('Incompatible size of sigma parametrisation!')
            end
        elseif(and(size(sigma_C,1)==1,size(sigma_C,2)==1))
            Sigma_C = repmat(sigma_C,size(Cdata));
        else
            error('Incompatible size of sigma parametrisation!')
        end
        
        SigmaStruct.Sigma_C = Sigma_C;
    end

end



function Sim = processSimulation(Sim, t_ind, thisData, type)

    % Mean Simulation
    Sim.m = Sim.mFine(t_ind,:);
    if strcmp(type, 'SCSH')
        Sim.C = Sim.CFine(t_ind,:);
        if isfield(Sim, 'Y_true')
            Sim.Y_true = Sim.Y_true(t_ind,:,:);
        end
    end
    
    if (size(thisData.condition,1) == 1)
        % No dose reponse experiment
        if (size(thisData.(type).time,1) ~= size(Sim.m,1))
            k = 0;
            oldT = nan;
            tmp_my = nan(size(thisData.(type).time,1), size(Sim.m,2));
            if strcmp(type, 'SCSH')
                tmp_Cy = nan(size(thisData.(type).time,1), size(Sim.C,2));
            end
            for j = 1 : size(thisData.(type).time,1)
                if (thisData.(type).time(j) ~= oldT)
                    k = k + 1; 
                end
                tmp_my(j,:) = Sim.m(k,:);
                if strcmp(type, 'SCSH')
                    tmp_Cy(j,:) = Sim.C(k,:);
                end
                oldT = thisData.(type).time(j);
            end
            Sim.m = tmp_my;
            if strcmp(type, 'SCSH')
                Sim.C = tmp_Cy;
            end
        end
    else
        % Dose response experiment
        tmp_my = nan(size(thisData.condition,1), size(Sim.m,2));
        thisUniqueCondition = unique(thisData.condition, 'rows');
        for j = 1 : size(thisData.condition,1) % number of conditions
            for iDose = 1 : size(thisUniqueCondition,1)
                if all(thisUniqueCondition(iDose,:) == thisData.condition(j,:))
                    tmp_my(j,:) =  Sim.m(iDose,:);
                end
            end
        end
        Sim.m = tmp_my;
    end

end