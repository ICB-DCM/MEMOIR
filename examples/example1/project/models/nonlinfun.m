function [status,t,x,y,sx,sy,s2x,s2y] = nonlinfun(tout,phi,kappa)
    status = 1;
    if(size(tout,2)>size(tout,1))
        tout = transpose(tout);
    end
    t = tout;
    x(:,1) = phi(1)*(1-tout) + phi(2)*tout;
    x(:,2) = phi(2)*(1-tout) + phi(1)*tout;
    
    y(:,1) = x(:,1);
    y(:,2) = x(:,2);
    
    sx(:,1,1) = (1-tout);
    sx(:,1,2) = tout;
    sx(:,2,1) = tout;
    sx(:,2,2) = (1-tout);
    
    sy(:,1,:) = sx(:,1,:);
    sy(:,2,:) = sx(:,2,:);
    
    s2x(:,1,1,1) = zeros(size(tout));
    s2x(:,1,1,2) = zeros(size(tout));
    s2x(:,1,2,1) = zeros(size(tout));
    s2x(:,1,2,2) = zeros(size(tout));
    
    s2x(:,2,1,1) = zeros(size(tout));
    s2x(:,2,1,2) = zeros(size(tout));
    s2x(:,2,2,1) = zeros(size(tout));
    s2x(:,2,2,2) = zeros(size(tout));
    
    s2y(:,1,:,:) = s2x(:,1,:,:);
    s2y(:,2,:,:) = s2x(:,2,:,:);
    
end