%function logL = logL_MEM_wo_grad(theta,t,D,options)
function varargout = logL_MEM_w_grad(varargin)

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
d = 10.^theta(5:end);

%% Evaluation of likelihood function
logL = 0;
dlogLdtheta = zeros(length(theta),1);
for i = 1:size(D,2)
    % Simulation of model
    y = x0*exp(-d(i)*t);
    dydtheta = zeros(length(t),length(theta));
    dydtheta(:,1) = exp(-d(i)*t)*(x0*log(10));
    dydtheta(:,4+i) = -x0*t.*exp(-d(i)*t)*(d(i)*log(10));

    % Log-likelihood
    logL = logL - 1/2*sum(log(2*pi*sigma_noise^2) + ((log(D(:,i)) - log(y))/sigma_noise).^2); % likelihood of data given single-cell parameter
    logL = logL - 1/2*(log(2*pi*sigma_d^2) + ((log(d(i)) - mu_d)/sigma_d).^2);     % likelihood of single cell parameters
    % Gradient of log-likelihood
    dlogLdtheta = dlogLdtheta + sum(bsxfun(@times,dydtheta,(log(D(:,i)) - log(y))./(y*sigma_noise^2)),1)';
    dlogLdtheta(4) = dlogLdtheta(4) - sum(1/sigma_noise - (log(D(:,i)) - log(y)).^2/sigma_noise^3)*(sigma_noise*log(10));
    dlogLdtheta(2) = dlogLdtheta(2) + (log(d(i)) - mu_d)/sigma_d^2;
    dlogLdtheta(3) = dlogLdtheta(3) - (1/sigma_d - (log(d(i)) - mu_d)^2/sigma_d^3)*(sigma_d*log(10));
    dlogLdtheta(4+i) = dlogLdtheta(4+i) - (log(d(i)) - mu_d)/(d(i)*sigma_d^2)*(d(i)*log(10));

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

