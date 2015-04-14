function [status,t,x,y,sx,sy] = simulate_RRE_mRNA_transfection(tout,phi,kappa)

% Transformation of parameters
theta = exp(phi);
    
persistent old_theta old_tout
persistent ncount
persistent old_status old_t old_x old_y old_sx old_sy

try
        if(and(all(tout==old_tout),all(theta==old_theta)))
                        status = old_status;
        t = old_t;
        x = old_x;
        y = old_y;
        sx = old_sx;
        sy = old_sy;
        ncount = ncount + 1;
    if ncount > 5
                status = -1;
    end
        return
    end
    
catch
end

ncount = 0;
options_cvode.cvodes_atol = 1e-08;
options_cvode.cvodes_rtol = 1e-08;
options_cvode.cvodes_maxsteps = 10000;

options_cvode.xs = ones(2,1);
options_cvode.ys = ones(1,1);
options_cvode.ps = ones(1,4);
options_cvode.tstart = theta(1);

options_cvode.qPositiveX = zeros(size(options_cvode.xs));
options_cvode.uNum = zeros(1,0);
options_cvode.vNum = zeros(1,2);
options_cvode.dvdxNum = zeros(2,2);
options_cvode.dvduNum = zeros(2,0);
options_cvode.dvdpNum = zeros(2,4);
options_cvode.suNum = zeros(0,4);
options_cvode.svNum = zeros(1,2);

% Simulation
[status,t,x,y,sx,sy] = RRE_mRNA_transfection(tout,theta(1:4),exp(kappa(1:0)),options_cvode);
x(tout<theta(1),1) = zeros(sum(tout<theta(1)),1);
sx(:,:,1) = -x*[-theta(3)*theta(4),0;theta(2),-theta(3)]';
sx(:,1,2) = zeros(length(tout),1);

% Transformation of sensitivities
sx = bsxfun(@times,sx,permute(theta(1:4),[3,2,1]));
sy = sx(:,2,:);


old_status = status;
old_t = t;
old_x = x;
old_y = y;
old_sx = sx;
old_sy = sy;
old_theta = theta;
old_tout = tout;
end
