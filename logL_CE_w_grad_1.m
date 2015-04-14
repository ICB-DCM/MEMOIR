%% CURRENTLY DEPRECATED
% function [logL] = logL_CE_w_grad_1(theta,Data,Model,options)
function varargout = logL_CE_w_grad_1(varargin)

%% Load
persistent tau

if isempty(tau)
 tau = clock;
end

%% Initialization
theta = varargin{1};
Data = varargin{2};
Model = varargin{3};

options.sign = 'positive';
options.cvo_type = 'diag-matrix-logarithm';
options.grad_ind = [1:length(theta)]';
options.tau_update = 3;
options.plot = 0;
if nargin == 4
    options = setdefault(varargin{4},options);
end

if etime(clock,tau) > options.tau_update
    options.plot = 1;
    tau = clock;
end

%% Evaluation of likelihood function
% Initialization
logL = 0;
dlogLdtheta = zeros(length(theta),1);

% Data types
% - Single-Cell Time Lapse (SCTL)
% - Population Average (PA)
% - Population Snapshot (PS)
data_type = {'SCTL','PA','PS'};

% Assignment of global variables
A = Model.A;
B = Model.B;
n_beta = Model.n_beta;
n_D = Model.n_D;
n_b = Model.n_b;
ind_beta = Model.ind_beta;
ind_D = Model.ind_D;
type_D = Model.type_D;

% Construct fixed effects and Covariance matrix
beta = theta(ind_beta);
[D,invD,dD,dinvD,HD,HinvD] = xi2D(theta(ind_D),type_D);

% Loop: Experimental conditions
for s = 1:length(Data)
    
    %% Construction of time vector
    t_s = [];
    for dtype = 1:length(data_type)
        if isfield(Data{s},data_type{dtype}) 
            t_s = union(eval(['Data{s}.' data_type{dtype} '.time']),t_s);
        end
    end
    
    %% Single cell time-lapse data
    if isfield(Data{s},'SCTL')
        % Evaluation of time index set
        [~,ind_t] = ismember(Data{s}.SCTL.time,t_s);
        
        % Empty plot
        if options.plot, figure(1);
            subplot(length(Data),2,2*(s-1)+1); hold off;
            subplot(length(Data),2,2*(s-1)+2); hold off;
        end
        
        % Loop: Indiviudal cells
        for i = 1:size(Data{s}.SCTL.Y,3)
            % Load single-cell data
            Ym_si = Data{s}.SCTL.Y(:,:,i);
            Sigma_si = Data{s}.SCTL.Sigma(:,:,i);
            ind = find(~isnan(Ym_si));

            % Construct single-cell parameter
            ind_b_si = Model.ind_b_si(s,i);
            b_si = theta(ind_b_si);
            phi_si = A*beta + B*b_si;
                        
            % Simulate model
            [~,~,~,y,~,sy] = Model.simulation{s}(t_s,phi_si);
            
            % Evaluation of likelihood
            Y_si = y(ind_t,:);
            logL = logL ...
                    - 0.5*sum(((Ym_si(ind) - Y_si(ind))./Sigma_si(ind)).^2) ...
                    - 0.5*log(det(D)) - 0.5*b_si'*invD*b_si;
            
            if nargout >= 2
                % beta
                for k = 1:length(ind_beta)
                    sY_si_k = sy(ind_t,:,:)*(A(:,k).*exp(phi_si));
                    dlogLdtheta(ind_beta(k)) = dlogLdtheta(ind_beta(k)) ...
                        + sum(((Ym_si(ind) - Y_si(ind))./Sigma_si(ind).^2).*sY_si_k(ind));
                end
                % b
                for k = 1:length(ind_b_si)
                    sY_si_k = sy(ind_t,:,:)*(B(:,k).*exp(phi_si));
                    dlogLdtheta(ind_b_si(k)) = dlogLdtheta(ind_b_si(k)) ...
                        + sum(((Ym_si(ind) - Y_si(ind))./Sigma_si(ind).^2).*sY_si_k(ind)) ...
                        - invD(k,:)*b_si;
                end
                % D
                for k = 1:length(ind_D)
                    dlogLdtheta(ind_D(k)) = dlogLdtheta(ind_D(k)) ...
                        - 0.5*trace(invD*dD(:,:,k)) ...
                        - 0.5*b_si'*dinvD(:,:,k)*b_si;
                end                
            end
            
            % Visulization
            if options.plot 
                subplot(length(Data),2,2*(s-1)+1);
                plot(t_s,y,'r'); hold on;
                plot(Data{s}.SCTL.time,Data{s}.SCTL.Y(:,:,i),'b-');
                
                subplot(length(Data),2,2*(s-1)+2);
                plot(Data{s}.SCTL.time,Data{s}.SCTL.Y(:,:,i)-y,'b-');hold on;
            end
        end
        % Visulization
        if options.plot 
            subplot(length(Data),2,2*(s-1)+1);
            xlabel('time'); ylabel('y and y^m');
            title(Data{s}.name)

            subplot(length(Data),2,2*(s-1)+2);
            xlabel('time'); ylabel('y^m-y');
            title(Data{s}.name)

            drawnow;
        end
    end
    
    % Population average data
    
    % Population snapshot data
    
end

%% Output
if nargout <= 1
    % One output
    switch  options.sign
        case 'positive'
            varargout{1} =  logL;
        case 'negative'
            varargout{1} = -logL;
    end
else
    % Two outputs
    switch  options.sign
        case 'positive'
            varargout{1} =  logL;
            varargout{2} =  dlogLdtheta(options.grad_ind);
        case 'negative'
            varargout{1} = -logL;
            varargout{2} = -dlogLdtheta(options.grad_ind);
    end
end



