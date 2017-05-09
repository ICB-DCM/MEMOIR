function [SP,my,Cy,dmydxi,dCydxi]  = getSimulationSCSH(xi,Model,Data,s)
    % GETSIMULATIONSCSH Summary of this function goes here
    %   Detailed explanation goes here
    % Simulation using sigma points
    
    % Set options for sigma point routine
    nderiv = 0.5 * (nargout-1) - 1;
    op_SP.nderiv = nderiv;
    op_SP.req = [1,1,0,0,0,1,0];
    op_SP.type_D = Model.type_D;
   
    %% Simulate with a loop over different doses
    % Initialize
    my = [];
    Cy = [];
    dmydxi = [];
    dCydxi = [];
    
    % Loop over doses
    thisUniqueCondition = unique(Data{s}.condition,'rows');
    for iDose = 1:size(thisUniqueCondition,1) % number of conditions
        % Simulate
        SP = getSigmaPointApp(...
            @(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCSH.time,phi,Data{s}.condition), ... = nonfun (in getSigmaPointApp)
            xi, ...
            Model.exp{s}, ... = estruct (in getSigmaPointApp)
            op_SP);
        
% Old code part with logarithm
%         tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))),1:size(SP.Cy,1),'UniformOutput',false);
%         [my,exp(SP.my+transpose([tmp{:}])/2)];
%         [my,exp(SP.my+transpose([tmp{:}])/2)];
        
        % Store means and variances
        my = [my; SP.my]; 
        Cy = [my; SP.Cy]; 
        
        % Store gradients of means and variances
        if(nderiv>0)
            dmydxi = [dmydxi; SP.dmydxi];
            dCydxi = [dCydxi; SP.dCydxi];
            
% Old code part, if logarithm was taken in simulateForSP:
%             nt = size(SP.dCydxi,1);
%             np = size(SP.dCydxi,4);
%             ny = size(SP.dCydxi,2);
%             dtmpdxi = arrayfun(@(x,y) diag(squeeze(SP.dCydxi(x,:,:,y))),repmat(1:nt,[np,1]),...
%                 repmat(transpose(1:np),[1,nt]),'UniformOutput',false);
%             dmydxi = bsxfun(@times,my,SP.dmydxi) ...
%                 + bsxfun(@times,my,permute(reshape([dtmpdxi{:}]/2,...
%                 [ny,np,nt]),[3,1,2]));
%             dmydxi(isnan(dmydxi)) = 0;

        else
            dmydxi = [dmydxi; zeros([size(my,1),size(my,2),length(xi)])];
            dCydxi = [dCydxi; zeros([size(my,1),size(my,2),length(xi)])];
        end
    end
    
    
    %% Post-process and clean-up
    
    % Kill the NANs, althoug nansum ist used later... (necessary?)
    my(isnan(my)) = 0;
    Cy(isnan(Cy)) = 0;
    if (nderiv>0)
        dmydxi(isnan(dmydxi)) = 0;
        dCydxi(isnan(dCydxi)) = 0;
    end
    
% Old code part, no clue how useful:
%     nt = size(SP.dCydxi,1);
%     np = size(SP.dCydxi,4);
%     ny = size(SP.dCydxi,2);
%     dtmpdxi = arrayfun(@(x,y) diag(squeeze(SP.dCydxi(x,:,:,y))),repmat(1:nt,[np,1]),...
%         repmat(transpose(1:np),[1,nt]),'UniformOutput',false);
%     dmydxi = bsxfun(@times,my,SP.dmydxi) ...
%         + bsxfun(@times,my,permute(reshape([dtmpdxi{:}]/2,...
%         [ny,np,nt]),[3,1,2]));
%     dmydxi(isnan(dmydxi)) = 0;
    
    
    % Post-processing of single cell snap-shot data
    if isfield(Model.exp{s},'SCSH_post_processing')
        if(nderiv==1)
            SP.dmydxi = zeros([size(SP.my) size(xi,1)]);
        end
        [my, Cy, dmydxi, dCydxi] = Model.exp{s}.SCSH_post_processing(my, Cy, dmydxi, dCydxi);
    end
    
end

