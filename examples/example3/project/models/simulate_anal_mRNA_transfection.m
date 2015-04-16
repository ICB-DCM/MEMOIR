%function [status,t,X,Y,SX,SY,S2X,S2Y] = simulate_anal_mRNA_transfection(t,phi,kappa)
function sol = simulate_anal_mRNA_transfection(t,phi,kappa,options)
    
if(isempty(options))
    options.sensi = 0;
end

    theta = exp(phi);
    
    t0 = min(theta(1),max(t));
    kTL_m0 = theta(2);
    beta = min(theta(3),4);
    delta = min(theta(4),4);
    
    t = t(:);
    
    % compute states
    X = [        exp(-delta*(t-t0)).*(t>t0),...
        kTL_m0*(exp( -beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
    % compute output
    Y = X(:,2);
    
    status = 0;
    
    sol.status = status;
    sol.x = X;
    sol.y = Y;
    
    % first order sensitivities required
    if options.sensi >= 1
        % state sensitivities
        sX(:,:,1) = -X*[-delta,0;kTL_m0,-beta]';
        sX(:,:,2) = [ zeros(size(t)),...
            (exp(-beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
        sX(:,:,3) = [ zeros(size(t)),...
            kTL_m0*((-(t-t0).*exp( -beta*(t-t0)))/(delta-beta) + (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        sX(:,:,4) = [-exp(-delta*(t-t0)).*(t-t0).*(t>t0),...
            kTL_m0*((+(t-t0).*exp(-delta*(t-t0)))/(delta-beta) - (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        SX = bsxfun(@times,sX,permute(theta,[3,2,1]));
        % output sensitivities
        SY = SX(:,2,:);
        sol.sx = SX;
        sol.sy = SY;
    end
    % second order sensitivities required
    if options.sensi >=2
        % state sensitivities
        s2X(:,:,1,1) = -sX(:,:,1)*[-delta,0;kTL_m0,-beta]'; % ok
        s2X(:,:,1,2) = -sX(:,:,2)*[-delta,0;kTL_m0,-beta]'-X*[0,0;1,0]'; %ok
        s2X(:,:,1,3) = -sX(:,:,3)*[-delta,0;kTL_m0,-beta]'-X*[0,0;0,-1]'; %ok
        s2X(:,:,1,4) = -sX(:,:,4)*[-delta,0;kTL_m0,-beta]'-X*[-1,0;0,0]'; %ok
        
        s2X(:,:,2,1) = -sX(:,:,2)*[-delta,0;kTL_m0,-beta]'-X*[0,0;1,0]'; %ok
        s2X(:,:,2,2) = zeros(size(X)); %ok
        s2X(:,:,2,3) = [ zeros(size(t)),...% ok
            ((-(t-t0).*exp( -beta*(t-t0)))/(delta-beta) + (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        s2X(:,:,2,4) = [ zeros(size(t)),...% ok
            ((+(t-t0).*exp(-delta*(t-t0)))/(delta-beta) - (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        
        s2X(:,:,3,1) = -sX(:,:,3)*[-delta,0;kTL_m0,-beta]'-X*[0,0;0,-1]'; %ok
        s2X(:,:,3,2) = [ zeros(size(t)),...% ok
            ((-(t-t0).*exp( -beta*(t-t0)))/(delta-beta) + (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        s2X(:,:,3,3) = [ zeros(size(t)),...
            kTL_m0*(((t-t0).^2.*exp( -beta*(t-t0)))/(delta-beta) ... %ok
            - 2*(t-t0).*exp(-beta*(t-t0))/((delta-beta)^2) ...
            + 2*(exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^3).*(t>t0)];
        s2X(:,:,3,4) = [ zeros(size(t)),... %ok
            kTL_m0*(((t-t0).*(exp(-beta*(t-t0)) + exp(-delta*(t-t0))))/((delta-beta)^2) ...
            - 2*(exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^3).*(t>t0)];
        
        s2X(:,:,4,1) = -sX(:,:,4)*[-delta,0;kTL_m0,-beta]'-X*[-1,0;0,0]'; %ok
        s2X(:,:,4,2) = [ zeros(size(t)),...% ok
            ((+(t-t0).*exp(-delta*(t-t0)))/(delta-beta) - (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        s2X(:,:,4,3) = [ zeros(size(t)),... %ok
            kTL_m0*(((t-t0).*(exp(-beta*(t-t0)) + exp(-delta*(t-t0))))/((delta-beta)^2) ...
            - 2*(exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^3).*(t>t0)];
        s2X(:,:,4,4) = [exp(-delta*(t-t0)).*(t-t0).^2.*(t>t0),...
            kTL_m0*(-((t-t0).^2.*exp(-delta*(t-t0)))/(delta-beta) ...
            - 2*(t-t0).*exp(-delta*(t-t0))/((delta-beta)^2) ...
            + 2*(exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^3).*(t>t0)];
        
        
        S2X = bsxfun(@times,s2X,permute(theta*theta',[4,3,2,1])) + bsxfun(@times,sX,permute(diag(theta),[4,3,2,1]));
        % output sensitivities
        S2Y = S2X(:,2,:,:);
        
        sol.s2x = S2X;
        sol.s2y = S2Y;
    end
    
    


    
    
    
    % options = odeset('reltol',1e-10,'abstol',1e-10);
    % [~,X] = ode15s(@(t,x) [-delta,0;kTL_m0,-beta]*x*(t.^10/(t0^10+t.^10)),t,[1;0],options);
    % X = X(:,2);
    
end
