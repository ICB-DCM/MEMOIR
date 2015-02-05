%clc;

%theta = [beta;logD];
theta = parameters.MS.MAP.par+0.1*randn(size(parameters.MS.MAP.par));

eps = 1e-5;
[logL,dlogLdtheta] = logL_CE_wo_grad_1(theta,Data,Model);

logL_i = zeros(length(theta),1);
for i = 1:length(theta)
    fprintf('.');
    theta_i = theta;
    theta_i(i) = theta_i(i) + eps;
    logL_i(i) = logL_CE_wo_grad_1(theta_i,Data,Model);
end
ddd = [(logL_i-logL)/eps,dlogLdtheta]
ddd(1:9,:)
