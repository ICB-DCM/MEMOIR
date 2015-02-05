%function logL = logL_MEM_m_w_grad(theta,t,D,options)
function varargout = logL_MEM_m_w_grad(varargin)

%% Initialization
theta = varargin{1};
t = varargin{2};
D = varargin{3};

options.sign = 'positive';
options.grad_ind = [1:length(theta)]';
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% Assignment of parameters
x0 = 10.^theta(1);
mu_d = theta(2);
sigma_d = 10.^theta(3);
sigma_noise = 10.^theta(4);
g = linspace(-5,5,5001);
log_d = mu_d + sigma_d*g;
d = exp(log_d);

%% Evaluation of likelihood function
% Simulation of model
y = x0*exp(-t*d);
dydtheta = zeros(length(t),length(g),length(theta));
dydtheta(:,:,1) = y*log(10);
dydtheta(:,:,2) = bsxfun(@times,y,-t*d);
%x0*exp(-t*exp(log_d))*(-t*dexp(log_d)/dmu_d)
%x0*exp(-t*exp(log_d))*(-t*dexp(mu_d + sigma_d*g)/dmu_d)
%x0*exp(-t*exp(log_d))*(-t* exp(mu_d + sigma_d*g))

dydtheta(:,:,3) = bsxfun(@times,y,bsxfun(@times,-t,d.*g))*(sigma_d*log(10));
%dydtheta(:,:,3) = bsxfun(@times,y,bsxfun(@times,-t,d.*g))*(sigma_d*log(10));
%x0*exp(-t*exp(log_d))*(-t*dexp(log_d)/dsigma_d)
%x0*exp(-t*exp(log_d))*(-t*dexp(mu_d + sigma_d*g)/dsigma_d)
%x0*exp(-t*exp(log_d))*(-t* exp(mu_d + sigma_d*g)*g)

% Evaluation of likelihood of log_d
%log_p_log_d = - 1/2*(log(2*pi*sigma_d^2) + ((log_d - mu_d)/sigma_d).^2);
%log_p_log_d = - 1/2*(log(2*pi*sigma_d^2) + ((mu_d + sigma_d*g - mu_d)/sigma_d).^2);
%log_p_log_d = - 1/2*(log(2*pi*sigma_d^2) + g.^2);
%log_p_log_d = - 1/2*(log(2*pi) + 2*log(sigma_d) + g.^2);
log_p_log_d = - 1/2*(log(2*pi*sigma_d^2) + g.^2);
% dlog_p_log_ddtheta = zeros(length(g),length(theta));
% dlog_p_log_ddtheta(:,3) = -log(10);

% Evaluation of likelihood
logL = 0;
dlogLdtheta = zeros(length(theta),1);
for i = 1:size(D,2)
    % Likelihood of the indivuel cells
    logL_i = - 1/2*sum(log(2*pi*sigma_noise^2) + (bsxfun(@minus,log(D(:,i)),log(y))/sigma_noise).^2,1);
    dlogL_idtheta = sum(bsxfun(@times,bsxfun(@minus,log(D(:,i)),log(y))./(y*sigma_noise^2),dydtheta),1);
    dlogL_idtheta(:,:,4) = dlogL_idtheta(:,:,4) ...
        - sum(1/sigma_noise - bsxfun(@minus,log(D(:,i)),log(y)).^2/sigma_noise^3,1)*(sigma_noise*log(10));
    dlogL_idtheta = squeeze(dlogL_idtheta);
    % Summation of log-likelihood
    f_i = exp(logL_i+log_p_log_d);
    df_idtheta = bsxfun(@times,dlogL_idtheta',exp(logL_i+log_p_log_d));
%    df_idtheta = bsxfun(@times,(dlogL_idtheta+dlog_p_log_ddtheta)',exp(logL_i+log_p_log_d));
    logL = logL + log(1/2*(f_i(1:end-1)+f_i(2:end))*diff(log_d)');
    dlogLdtheta = dlogLdtheta ...
        + ((df_idtheta(:,1:end-1)+df_idtheta(:,2:end))*diff(log_d)')...
                     ./((f_i(1:end-1)+f_i(2:end))*diff(log_d)');
end

%%
if max(isnan(dlogLdtheta))
    logL = -inf;
    dlogLdtheta = zeros(length(theta),1);
end

%% Output
switch  options.sign
    case 'positive'
        varargout{1} =  logL;
        varargout{2} =  dlogLdtheta(options.grad_ind);
    case 'negative'
        varargout{1} = -logL;
        varargout{2} = -dlogLdtheta(options.grad_ind);
end

