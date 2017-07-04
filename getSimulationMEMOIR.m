function getSimulationMEMOIR(xi, Model, Data, ind_exp, ind_xi, options)
% This routine plots a number of specified experiments (ind_exp) with a 
% fine time grid and analyzes relative sensitivities for a given parameter 
% set (ind_xi)

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
        
        doseResponse = (size(unique(Data{s}.condition,'rows'),1) ~= 1);
        
        if doseResponse
            % If the current experiment is a dose-response experiment
            
        else
        % --- no dose response experiment ---------------------------------
            %% Simulation of the model
            % Set options for sigma point routine
            op_SP.nderiv = options.sensi;
            op_SP.req = req;
            op_SP.type_D = Model.type_D;
            op_SP.approx = 'sp';

            % Call simulation
            SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model, t_sim, phi, Data{s}.condition, Model.exp{s}.scale), ...  = nonfun (in getSigmaPointApp)
                xi, ...
                Model.exp{s}, ... = estruct (in getSigmaPointApp)
                op_SP);
            
            
            %% Post-processing
            % Store the simulation results and apply scaling
            switch Model.exp{s}.scale
                case 'log'
                    % TBD!
                    tmp = arrayfun(@(x) diag(squeeze(SP.my(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                    my = exp(SP.my + transpose([tmp{:}])/2);
                    if strcmp(exp_type, 'SCSH')
                        tmpCy = bsxfun(@plus, repmat(tmp, 1, 1, size(SP.Cy,3)), permute(repmat(tmp, 1, 1, size(SP.Cy,3)), [1,3,2]));
                        Cy = exp(tmpCy) .*  (exp(SP.Cy) - ones(size(SP.Cy)));
                    end

                case 'log10'
                    % TBD!
                    tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                    my = 10.^(SP.my + transpose([tmp{:}])/2);
                    if strcmp(exp_type, 'SCSH')
                        tmpCy = bsxfun(@plus, repmat(tmp, 1, 1, size(SP.Cy,3)), permute(repmat(tmp, 1, 1, size(SP.Cy,3)), [1,3,2]));
                        Cy = 10.^(tmpCy) .*  (10.^(SP.Cy) - ones(size(SP.Cy)));
                    end

                case 'lin'
                    my = SP.my;
                    if strcmp(exp_type, 'SCSH')
                        Cy = SP.Cy;
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

                    case 'log10'
                        % To be checked!
                        tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                        my = 10.^(SP.my + transpose([tmp{:}])/2);

                    case 'lin'
                        dmydxi = SP.dmydxi;
                        % dCydxi = SP.dCydxi;
                end
            else
                dmydxi = zeros([size(my,1),size(my,2),length(xi)]);
                % dCydxi = zeros([size(my,1),size(my,2),size(my,2),length(xi)]);
            end
            
            % Apply additional user-defined post-processing
            if isfield(Model.exp{s}, post_proc)
                if(options.sensi > 0)
                    SP.dmydxi = zeros([size(SP.my) size(xi,1)]);
                end
                [my,dmydxi] = feval(Model.exp{s}.(post_proc), my, dmydxi);
            end

            
            %% Plotting of fine simulation
            % Assign values for plotting
            Sim.mFine = my;
            Sim.t = unique(t_sim);
            SigmaStruct = processSigma(Data{s}, Sim.mFine, [], exp_type);
            Sim.Sigma_m = SigmaStruct.Sigma_m;
            if strcmp(exp_type, 'SCSH')
                Sim.CFine = Cy;
                SigmaStruct = processSigma(Data{s}, Sim.mFine, Cy, exp_type);
                Sim.Sigma_C = SigmaStruct.Sigma_C;
            end
            
            % Process values in Sim.m for residual plot
            Sim = processSimulation(Sim, t_ind, Data{s}, exp_type);
            
            % Plotting
            Model.exp{s}.plot(Data{s}, Sim, s);
            
        % --- End of no dose response experiment --------------------------
        end
        
        % Clean up
        clear Sim;
    % --- End of loop over experiments ------------------------------------
    end
    
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
        
    if strcmp(type, 'SCSH')
        % Single-Cell snapshot sigma
        sigma_C = thisData.(type).Sigma_C(1,:);
        
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
    end
    
    if (size(thisData.(type).time,1) ~= size(Sim.m,1))
        k = 0;
        oldT = nan;
        tmp_my = nan(size(thisData.(type).time,1), size(Sim.m,2));
        for j = 1 : size(thisData.(type).time,1)
            if (thisData.(type).time(j) ~= oldT)
                k = k + 1; 
            end
            tmp_my(j,:) = Sim.m(k,:);
            oldT = thisData.(type).time(j);
        end
        Sim.m = tmp_my;
    end

end