clc;

theta = parameters.MS.par(:,1)+0.1*rand(size(parameters.MS.par(:,1)));

eps = 1e-5;
[logL,dlogLdtheta] = logL_CE_w_grad_2(theta,Data,Model);

logL_i = zeros(length(theta),1);
for i = 1:length(theta)
    theta_i = theta;
    theta_i(i) = theta_i(i) + eps;
    logL_i(i) = logL_CE_w_grad_2(theta_i,Data,Model);
end

ddd = [(logL_i-logL)/eps,dlogLdtheta]
