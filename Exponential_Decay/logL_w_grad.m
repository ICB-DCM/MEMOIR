%function logL = logL_wo_grad(theta,t,D,options)
function varargout = logL_w_grad(varargin)

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
d = 10.^theta(2);
sigma_noise = 10.^theta(3);

%% Evaluation of likelihood function
% Simulation of model
y = x0*exp(-d*t);
dydtheta = [       exp(-d*t)*(x0*log(10)),...
            -x0*t.*exp(-d*t)*(d *log(10))];
% Log-likelihood
logL = - 1/2*sum(sum(log(2*pi*sigma_noise^2) + (bsxfun(@minus,log(D),log(y))/sigma_noise).^2));
dlogLdtheta(1,1) = sum(sum(bsxfun(@times,bsxfun(@minus,log(D),log(y)),dydtheta(:,1)./y)/sigma_noise^2));
dlogLdtheta(2,1) = sum(sum(bsxfun(@times,bsxfun(@minus,log(D),log(y)),dydtheta(:,2)./y)/sigma_noise^2));
dlogLdtheta(3,1) = sum(sum(- 1 + (bsxfun(@minus,log(D),log(y))/sigma_noise).^2))*log(10);

%% Output
switch  options.sign
    case 'positive'
        varargout{1} =  logL;
        varargout{2} =  dlogLdtheta(options.grad_ind);
    case 'negative'
        varargout{1} = -logL;
        varargout{2} = -dlogLdtheta(options.grad_ind);
end

