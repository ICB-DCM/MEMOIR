%function logL = logL_MEM_wo_grad(theta,t,D,options)
function varargout = logL_MEM_wo_grad(varargin)

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
d = 10.^theta(5:end);

%% Evaluation of likelihood function
logL = 0;
for i = 1:size(D,2)
    % Simulation of model
    y = x0*exp(-d(i)*t);
    % Summation of log-likelihood
    logL = logL - 1/2*sum(log(2*pi*sigma_noise^2) + ((log(D(:,i)) - log(y))/sigma_noise).^2); % likelihood of data given single-cell parameter
    logL = logL - 1/2*(log(2*pi*sigma_d^2) + ((log(d(i)) - mu_d)/sigma_d).^2);     % likelihood of single cell parameters
end

%% Output
switch  options.sign
    case 'positive'
        varargout{1} =  logL;
    case 'negative'
        varargout{1} = -logL;
end

