function sol = nonlinfun(tout,phi,kappa,options)
    sol.status = 1;
    if(size(tout,2)>size(tout,1))
        tout = transpose(tout);
    end
    sol.t = tout;
    
    % compute states
    x(:,1) = phi(1)*(1-tout) + phi(2)*tout;
    x(:,2) = phi(2)*(1-tout) + phi(1)*tout;
    
    % compute output
    y(:,1) = x(:,1);
    y(:,2) = x(:,2);
    
    sol.x = x;
    sol.y = y;
    
    % first order sensitivities required
    if(options.sensi >= 1)
        % state sensitivities
        sx(:,1,1) = (1-tout);
        sx(:,1,2) = tout;
        sx(:,2,1) = tout;
        sx(:,2,2) = (1-tout);
        % output sensitivities
        sy(:,1,:) = sx(:,1,:);
        sy(:,2,:) = sx(:,2,:);
        sol.sx = sx;
        sol.sy = sy; 
    end
    
    % second order sensitivities required
    if(options.sensi >= 2)
        % state sensitivities
        s2x(:,1,1,1) = zeros(size(tout));
        s2x(:,1,1,2) = zeros(size(tout));
        s2x(:,1,2,1) = zeros(size(tout));
        s2x(:,1,2,2) = zeros(size(tout));
        
        s2x(:,2,1,1) = zeros(size(tout));
        s2x(:,2,1,2) = zeros(size(tout));
        s2x(:,2,2,1) = zeros(size(tout));
        s2x(:,2,2,2) = zeros(size(tout));
       
        % output sensitivities
        s2y(:,1,:,:) = s2x(:,1,:,:);
        s2y(:,2,:,:) = s2x(:,2,:,:);

        sol.s2x = s2x;
        sol.s2y = s2y;
    end
end