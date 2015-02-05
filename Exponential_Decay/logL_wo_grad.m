%function logL = logL_wo_grad(theta,t,D,options)
function varargout = logL_wo_grad(varargin)

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
d = 10.^theta(2);
sigma_noise = 10.^theta(3);

%% Evaluation of likelihood function
% Simulation of model
y = x0*exp(-d*t);
% Summation of log-likelihood
logL = - 1/2*sum(sum(log(2*pi*sigma_noise^2) + (bsxfun(@minus,log(D),log(y))/sigma_noise).^2));

%% Output
switch  options.sign
    case 'positive'
        varargout{1} =  logL;
    case 'negative'
        varargout{1} = -logL;
end

