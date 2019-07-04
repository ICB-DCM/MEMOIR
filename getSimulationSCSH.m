function [SP,my,Cy,dmydxi,dCydxi]  = getSimulationSCSH(xi,Model,Data,s,options)
    % GETSIMULATIONSCSH Summary of this function goes here
    %   Detailed explanation goes here
    % Simulation using sigma points
    
    % Set options for sigma point routine
    nderiv = 0.5 * (nargout-1) - 1;
    op_SP.nderiv = nderiv;
    op_SP.req = [1,1,0,0,0,1,0];
    op_SP.type_D = Model.type_D;
    op_SP.approx = options.approx;
    if isfield(options, 'samples')
        op_SP.samples = options.samples;
    end
    
    %% Simulate with a loop over different doses
    % Initialize
    my = [];
    Cy = [];
    dmydxi = [];
    dCydxi = [];
    
    % Loop over doses
    thisUniqueCondition = unique(Data{s}.condition,'rows');
    for iDose = 1:size(thisUniqueCondition,1)
    % === Loop over doses =================================================
        % Simulate
        SP = getSigmaPointApp(...
            @(phi) simulateForSP(Model.exp{s}.model, Data{s}.SCSH.time, phi, thisUniqueCondition(iDose,:), Model.exp{s}.scale), ... = nonfun (in getSigmaPointApp)
            xi, ...
            Model.exp{s}, ... = estruct (in getSigmaPointApp)
            op_SP);
        
        % Store the simulation results
        switch Model.exp{s}.scale
            case 'log'
                tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))), 1:size(SP.Cy,1),'UniformOutput',false);
                tmp = transpose([tmp{:}])/2;
                my_new = exp(SP.my + tmp);
                my = [my; my_new];
                tmpCy = bsxfun(@plus, repmat(tmp, 1, 1, size(SP.Cy,3)), permute(repmat(tmp, 1, 1, size(SP.Cy,3)), [1,3,2]));
                Cy_new = exp(tmpCy) .*  (exp(SP.Cy) - ones(size(SP.Cy)));
                Cy = [Cy; Cy_new];
                
            case 'log10'
                tmp = arrayfun(@(x) diag(squeeze(SP.my(x,:,:))), 1:size(SP.my,1),'UniformOutput',false);
                tmp = transpose([tmp{:}])/2;
                my = [my; 10.^(SP.my + tmp)];
                tmpCy = bsxfun(@plus, repmat(tmp, 1, 1, size(SP.Cy,3)), permute(repmat(tmp, 1, 1, size(SP.Cy,3)), [1,3,2]));
                Cy = [Cy; 10.^(tmpCy) .*  (10.^(SP.Cy) - ones(size(SP.Cy)))];
                
            case 'lin'
                if strcmp(op_SP.approx, 'pa only')
                    my = [my; SP.my];
                else
                    my = [my; SP.my];
                    Cy = [Cy; SP.Cy];
                end
        end

        % Store gradients of means and variances
        if(nderiv>0)
            switch Model.exp{s}.scale
                case 'log'
                    nt = size(SP.dCydxi,1);
                    np = size(SP.dCydxi,4);
                    ny = size(SP.dCydxi,2);
                    dtmpdxi = arrayfun(@(x,y) diag(squeeze(SP.dCydxi(x,:,:,y))),repmat(1:nt,[np,1]),...
                        repmat(transpose(1:np),[1,nt]),'UniformOutput',false);
                    dtmpdxi = permute(reshape([dtmpdxi{:}]/2, [ny,np,nt]),[3,1,2]);
                    dmydxi = [dmydxi; bsxfun(@times,my_new,SP.dmydxi) ...
                        + bsxfun(@times,my_new, dtmpdxi)];
                    dmydxi(isnan(dmydxi)) = 0;
                    dtmpCydxi = bsxfun(@plus, repmat(permute(dtmpdxi, [1,2,4,3]),...
                        1, 1, size(SP.Cy,3), 1), repmat(permute(dtmpdxi, [1,4,2,3]),...
                        1, size(SP.Cy,3), 1, 1));
                    dCydxi = [dCydxi; bsxfun(@times, exp(tmpCy) .* (exp(SP.Cy) - ones(size(SP.Cy))), dtmpCydxi)...
                        + bsxfun(@times, exp(tmpCy), bsxfun(@times, exp(SP.Cy), SP.dCydxi))];

                case 'log10'
                    nt = size(SP.dCydxi,1);
                    np = size(SP.dCydxi,4);
                    ny = size(SP.dCydxi,2);
                    dtmpdxi = arrayfun(@(x,y) diag(squeeze(SP.dCydxi(x,:,:,y))),repmat(1:nt,[np,1]),...
                        repmat(transpose(1:np),[1,nt]),'UniformOutput',false);
                    dmydxi = bsxfun(@times,my,SP.dmydxi) ...
                        + bsxfun(@times,my,permute(reshape([dtmpdxi{:}]/2,...
                        [ny,np,nt]),[3,1,2]));
                    dmydxi(isnan(dmydxi)) = 0;

                case 'lin'
                    if strcmp(op_SP.approx, 'pa only')
                        dmydxi = [dmydxi; SP.dmydxi];
                    else
                        dmydxi = [dmydxi; SP.dmydxi];
                        dCydxi = [dCydxi; SP.dCydxi];
                    end
            end
        else
            dmydxi = [dmydxi; zeros([size(my,1),size(my,2),length(xi)])];
            dCydxi = [dCydxi; zeros([size(my,1),size(my,2),size(my,2),length(xi)])];
        end
    % === Loop over doses ends ============================================ 
    end
    
    
    %% Post-process and clean-up
    if strcmp(op_SP.approx, 'pa only')
        Cy = zeros(size(Data{s}.SCSH.C));
        dCydxi = zeros([size(Data{s}.SCSH.C, 1), size(Data{s}.SCSH.C, 2), length(xi)]);
        SP.Cy = Cy;
        SP.dCydxi = dCydxi;
    end
    
    % Kill the NANs, althoug nansum ist used later... (necessary?)
    my(isnan(my)) = 0;
    Cy(isnan(Cy)) = 0;
    if (nderiv>0)
        dmydxi(isnan(dmydxi)) = 0;
        dCydxi(isnan(dCydxi)) = 0;
    end
    
    % Post-processing of single cell snap-shot data
    if isfield(Model.exp{s},'SCSH_post_processing')
        if(nderiv==1)
            SP.dmydxi = zeros([size(SP.my) size(xi,1)]);
            SP.dCydxi = zeros([size(Cy) size(xi,1)]);
        end
        [my, Cy, dmydxi, dCydxi] = Model.exp{s}.SCSH_post_processing(my, Cy, dmydxi, dCydxi);
    end
    
end
