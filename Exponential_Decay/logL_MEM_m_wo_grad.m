%function logL = logL_MEM_wo_grad(theta,t,D,options)
function varargout = logL_MEM_m_wo_grad(varargin)

%% Initialization
theta = varargin{1};
t = varargin{2};
D = varargin{3};

options.sign = 'positive';
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% Assignment of parameters
x0 = 10.^theta(1);
mu_d = theta(2);
sigma_d = 10.^theta(3);
sigma_noise = 10.^theta(4);
g = linspace(-5,5,501);
log_d = mu_d + sigma_d*g;
d = exp(log_d);

%% Evaluation of likelihood function
% Simulation of model
y = x0*exp(-t*d);
% Evaluation of likelihood of log_d
%log_p_log_d = - 1/2*(log(2*pi*sigma_d^2) + ((log_d - mu_d)/sigma_d).^2);
log_p_log_d = - 1/2*(log(2*pi*sigma_d^2) + g.^2);

% Evaluation of likelihood
logL = 0;
for i = 1:size(D,2)
    logL_i = - 1/2*sum(log(2*pi*sigma_noise^2) + (bsxfun(@minus,log(D(:,i)),log(y))/sigma_noise).^2,1);
    % Summation of log-likelihood
    f_i = exp(logL_i+log_p_log_d);
    logL = logL + log(0.5*(f_i(1:end-1)+f_i(2:end))*diff(log_d)');
end

%% Output
switch  options.sign
    case 'positive'
        varargout{1} =  logL;
    case 'negative'
        varargout{1} = -logL;
end

