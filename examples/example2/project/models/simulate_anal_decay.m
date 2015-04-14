%function [status,t,X,Y,SX,SY,S2X,S2Y] = simulate_anal_decay(t,phi,kappa)
function sol = simulate_anal_decay(t,phi,kappa,options)
    
if(isempty(options))
    options.sensi = 0;
end

    theta = exp(phi);

    m0 = theta(1);
    delta  = theta(2);
    
    t = t(:);

    % Simulation
    X(:,1) = m0*exp(-delta*t);
    Y = X;
    if options.sensi >= 1
        sX = zeros(numel(t),1,2);
        sX(:,:,1) = exp(-delta*t);
        sX(:,:,2) = -t.*m0.*exp(-delta*t);

        SX = bsxfun(@times,sX,permute(theta,[3,2,1]));

        SY = SX;
    end
    if options.sensi >=2
        s2X = zeros(numel(t),1,2,2);
        s2X(:,:,1,1) = 0;
        s2X(:,:,1,2) = -t.*exp(-delta*t);
        s2X(:,:,2,1) = -t.*exp(-delta*t);
        s2X(:,:,2,2) = t.*t.*m0.*exp(-delta*t);

        S2X = bsxfun(@times,s2X,permute(theta*theta',[4,3,2,1])) + bsxfun(@times,sX,permute(diag(theta),[4,3,2,1]));
        
        S2Y = S2X;
    end
    
    status = 0;
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
end
