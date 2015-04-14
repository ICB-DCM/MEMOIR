function sol = nonlinfun(tout,phi,kappa,options)
    sol.status = 1;
    if(size(tout,2)>size(tout,1))
        tout = transpose(tout);
    end
    sol.t = tout;
    sol.x(:,1) = phi(1)*(1-tout) + phi(2)*tout;
    sol.x(:,2) = phi(2)*(1-tout) + phi(1)*tout;
    
    sol.y(:,1) = sol.x(:,1);
    sol.y(:,2) = sol.x(:,2);
    sol.root = double.empty([0,1]);
    sol.rootval = double.empty([0,1]);
    
    if(options.sensi>0)
        sol.sx(:,1,1) = (1-tout);
        sol.sx(:,1,2) = tout;
        sol.sx(:,2,1) = tout;
        sol.sx(:,2,2) = (1-tout);
        
        
        sol.sy(:,1,:) = sol.sx(:,1,:);
        sol.sy(:,2,:) = sol.sx(:,2,:);
        
        sol.sroot = double.empty([0,1,2]);
        sol.srootval = double.empty([0,1,2]);
    end
    if(options.sensi>1)
        sol.s2x(:,1,1,1) = zeros(size(tout));
        sol.s2x(:,1,1,2) = zeros(size(tout));
        sol.s2x(:,1,2,1) = zeros(size(tout));
        sol.s2x(:,1,2,2) = zeros(size(tout));
        
        sol.s2x(:,2,1,1) = zeros(size(tout));
        sol.s2x(:,2,1,2) = zeros(size(tout));
        sol.s2x(:,2,2,1) = zeros(size(tout));
        sol.s2x(:,2,2,2) = zeros(size(tout));
        
        sol.s2y(:,1,:,:) = sol.s2x(:,1,:,:);
        sol.s2y(:,2,:,:) = sol.s2x(:,2,:,:);
        sol.s2root = double.empty([0,1,2,2]);
        sol.s2rootval = double.empty([0,1,2,2]);
    end
    
end