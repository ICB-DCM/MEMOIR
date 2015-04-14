function [status,t,x,y,sx,sy,s2x,s2y] = simulate_RRE_mRNA_transfection_deg(tout,phi,kappa)

tout = tout(:);

% Transformation of parameters
theta = exp(phi);

% Simulation
x = theta(2)*exp(-theta(1)*tout);
y = x;


sx(:,1,1) = -theta(2)*tout.*exp(-theta(1)*tout);
sx(:,1,2) = exp(-theta(1)*tout);

s2x(:,1,1,1) = theta(2)*tout.^2.*exp(-theta(1)*tout);
s2x(:,1,1,2) = -tout.*exp(-theta(1)*tout);
s2x(:,1,2,1) = -tout.*exp(-theta(1)*tout);
s2x(:,1,2,2) = zeros(x);

% Transformation of sensitivities
sx = bsxfun(@times,sx,permute(theta,[3,2,1]));
sy = sx;

s2x = bsxfun(@times,s2x,permute(theta*theta',[4,3,2,1])) + bsxfun(@times,sx,permute(diag([1,1]),[4,3,2,1])); 

%
t = tout;

status = [];

end
