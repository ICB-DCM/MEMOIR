function [SP,my,Cy,dmydxi,dCydxi]  = getSimulationSCSH(xi,Model,Data,s )
    %GETSIMULATIONPA Summary of this function goes here
    %   Detailed explanation goes here
    % Simulation using sigma points
    
    % Simulation using sigma points
    nderiv = nargout-1;
    op_SP.nderiv = nderiv;
    op_SP.req = [1,1,0,0,0,1,0];
    op_SP.type_D = Model.type_D;
    SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCSH.time,phi,Data{s}.condition),xi,Model.exp{s},op_SP);
    %     [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(x)getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.PA.time,phi,Data{s}.condition),x,Model.exp{s},op_SP),1e-4,'my','dmydxi')
    
    % we log transformed in simulateForSP, now we need to compute the mean
    % from a lognormal density
    tmp = arrayfun(@(x) diag(squeeze(SP.Cy(x,:,:))),1:size(SP.Cy,1),'UniformOutput',false);
    my = exp(SP.my+transpose([tmp{:}])/2);
    my(isnan(my)) = 0;
    
    nt = size(SP.dCydxi,1);
    np = size(SP.dCydxi,4);
    ny = size(SP.dCydxi,2);
    dtmpdxi = arrayfun(@(x,y) diag(squeeze(SP.dCydxi(x,:,:,y))),repmat(1:nt,[np,1]),...
        repmat(transpose(1:np),[1,nt]),'UniformOutput',false);
    dmydxi = bsxfun(@times,my,SP.dmydxi) ...
        + bsxfun(@times,my,permute(reshape([dtmpdxi{:}]/2,...
        [ny,np,nt]),[3,1,2]));
    dmydxi(isnan(dmydxi)) = 0;
    
    
end
