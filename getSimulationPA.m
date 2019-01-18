function [SP,my,dmydxi]  = getSimulationPA(xi,Model,Data,s,options)
    %GETSIMULATIONPA Summary of this function goes here
    %   Detailed explanation goes here
    % Simulation using sigma points
    
    % Set options for sigma point routine
    nderiv = nargout-2;
    op_SP.nderiv = nderiv;
    op_SP.req = [1,0,0,0,0,1,0]; % [1,1,0,0,0,1,0];
    op_SP.type_D = Model.type_D;
    
    op_SP.approx = options.approx;
    if isfield(options, 'samples')
        op_SP.samples = options.samples;
    end

    %% Simulate with a loop over different doses
    % Initialize
    my = [];
    dmydxi = [];
    
    % Loop over doses
    thisUniqueCondition = unique(Data{s}.condition,'rows');
    for iDose = 1:size(thisUniqueCondition,1)
    % === Loop over doses =================================================
        % Simulate
        SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model, Data{s}.PA.time, phi, thisUniqueCondition(iDose,:), Model.exp{s}.scale), ...  = nonfun (in getSigmaPointApp)
            xi, ...
            Model.exp{s}, ... = estruct (in getSigmaPointApp)
            op_SP);

        % Store the simulation results
        switch Model.exp{s}.scale
            case 'log'
                % TBD!
                tmp = arrayfun(@(x) diag(squeeze(SP.my(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                my = [my; exp(SP.my + transpose([tmp{:}])/2)];

            case 'log10'
                % TBD!
                tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                my = [my; 10.^(SP.my + transpose([tmp{:}])/2)];
                
            case 'lin'
                my = [my; SP.my];
        end
        
        % Store gradients of means and variances
        if(nderiv>0)
            switch Model.exp{s}.scale
                case 'log'
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
                    my = [my, 10.^(SP.my + transpose([tmp{:}])/2)];

                case 'lin'
                    dmydxi = [dmydxi; SP.dmydxi]; 
            end
        else
            dmydxi = [dmydxi;zeros([size(my,1),size(my,2),length(xi)])];
        end
    % === Loop over doses ends ============================================        
    end

        
    %% Post-process and clean-up
    
    % Kill the NANs, althoug nansum ist used later... (necessary?)
    my(isnan(my)) = 0;
    if (nderiv>0)
        dmydxi(isnan(dmydxi)) = 0;
    end
    
    
    % Post-processing of population average data
    if isfield(Model.exp{s},'PA_post_processing')
        if(nderiv==1)
            SP.dmydxi = zeros([size(SP.my) size(xi,1)]);
        end
        [my,dmydxi] = Model.exp{s}.PA_post_processing(my, dmydxi, xi);
    end
    if isfield(Model.exp{s},'PA_post_processing_SP')
        SP = Model.exp{s}.PA_post_processing_SP(SP);
    end
end
