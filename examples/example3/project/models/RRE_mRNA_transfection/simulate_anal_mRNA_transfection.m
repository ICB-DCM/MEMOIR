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
    
    
    % theta = (t0,kTL_m0,beta,delta)
    
    % dm/dt = - delta*m       => m = exp(-delta*t)*m0
    % dG/dt = kTL*m - beta*G
    %
    % syms G(t)
    % G_sol = dsolve(diff(G) == kTL*exp(-delta*t)*m0 - beta*G, G(0) == 0)
    % G_sol = simplify(G_sol)
    
    % Simulation
    X = [        exp(-delta*(t-t0)).*(t>t0),...
        kTL_m0*(exp( -beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
    Y = X(:,2);
    if options.sensi >= 1
        sX(:,:,1) = -X*[-delta,0;kTL_m0,-beta]';
        sX(:,:,2) = [ zeros(size(t)),...
            (exp(-beta*(t-t0)) - exp(-delta*(t-t0)))/(delta-beta).*(t>t0)];
        sX(:,:,3) = [ zeros(size(t)),...
            kTL_m0*((-(t-t0).*exp( -beta*(t-t0)))/(delta-beta) + (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        sX(:,:,4) = [-exp(-delta*(t-t0)).*(t-t0).*(t>t0),...
            kTL_m0*((+(t-t0).*exp(-delta*(t-t0)))/(delta-beta) - (exp(-beta*(t-t0))-exp(-delta*(t-t0)))/(delta-beta)^2).*(t>t0)];
        
        
        SX = bsxfun(@times,sX,permute(theta,[3,2,1]));
        
        SY = SX(:,2,:);
    end
    if options.sensi >=2
        
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
        
        S2Y = S2X(:,2,:,:);
    end
    
    
    status = 1;
    
    sol.status = status;
    sol.x = X;
    sol.y = Y;
    sol.root = zeros(0,1);
    sol.rootval = zeros(0,1);
    if options.sensi >= 1
        sol.sx = SX;
        sol.sy = SY;
        sol.sroot = zeros(0,1,length(theta));
        sol.srootval = zeros(0,1,length(theta));
    end
    if options.sensi >= 2
        sol.s2x = S2X;
        sol.s2y = S2Y;
        sol.s2root = zeros(0,1,length(theta),length(theta));
        sol.s2rootval = zeros(0,1,length(theta),length(theta));
    end

    
    
    
    % options = odeset('reltol',1e-10,'abstol',1e-10);
    % [~,X] = ode15s(@(t,x) [-delta,0;kTL_m0,-beta]*x*(t.^10/(t0^10+t.^10)),t,[1;0],options);
    % X = X(:,2);
    
end
