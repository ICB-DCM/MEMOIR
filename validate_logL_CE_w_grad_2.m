% logL_CE_w_grad_2 is the main likelihood function of the MEMToolbox and
% supports classical mixed effect models aswell as sigma-point based
% moment-approximation schemes to the solution of population balance
% equations. The likelihood function allows for the parallel estimation of
% parameters from multiple different experimental setups/conditions
%
% USAGE:
% ======
% P = logL_CE_w_grad_2(xi,Data,Model,options)
%
% INPUTS:
% =======
% xi ... parameter of for which the logLikelihood value is to be computed
%   the ordering must correspond to what is given in Model.sym.xi
% Data ... data for which the logLikelihood is computed this should be a cell
%   array (one cell per experiment) with some of the following fields
%   .SCTL ... Single-Cell-Time-Lapse; must have the following subfields
%       .time ... vector of time-points at which measurements were taken
%       .Y ... measurements, 3rd order tensor [time,observable,cells]
%       .T ... (optional) events, 3rd order tensor [time,event,cells]
%   .SCTLstat ... Single-Cell-Time-Lapse-Statistics; must have the following subfields
%       .time ... vector of time-points at which measurements were taken
%       .mz ... full state vector of all observables at all times
%       .Sigma_mz ... standard deviation of mz
%       .Cz ... full state covariance, this includes state covariance and
%       cross-correlation
%       .Sigma_Cz ... standard deviation of Cz
%   .SCSH ... Single-Cell-Snap-Shot; must have the following subfields
%       .time ... vector of time-points at which measurements were taken
%       .m ... mean of measurements, matrix [time,observable]
%       .Sigma_m ... standard deviation of mean, matrix [time,observable]
%       .C ... covariance of measurements, 3rd order tensor [time,observable,observable]
%       .Sigma_C ... standard deviation of covariance, 3rd order tensor [time,observable,observable]
%   .PA ... Population-Average; must have the following subfields
%       .time ... vector of time-points at which measurements were taken
%       .m ... mean of measurements, matrix [time,observable]
%       .Sigma_m ... standard deviation of mean, matrix [time,observable]
%   furthermore it must have the following fields
%   .condition ... vector of experiemental conditions which is passed to
%       the model file
%   .name ... string for the name of the experiment
%   .measurands ... cell array of labels for every experiment
% Model ... model definition ideally generated via make_model and
%   complete_model. must have the following fields
%   .type_D ... string specifying the parametrisation of the covariance
%       matrix for the random effects. either
%       'diag-matrix-logarithm' for diagonal matrix with log. paramet. or
%       'matrix-logarithm' for full matrix with log. paramet. or
%   .integration ... flag indicating whether integration for classical
%       mixed effect models via laplace approximation should be applied.
%       only applicable for SCTL data.
%   .penalty ... flag indicating whether additional penalty terms for
%       synchronisation of parameters across experiments should be applied.
%       only applicable for SCTL data.
%   .prior ... cell array containing prior information for the optimization
%       parameters. the length of the cell array should not exceed the
%       length of the optimization parameter. the prior is applied when
%       there exist fields .mu and .std. these two parameters define a
%       quadratic function centered around .mu multiplied with 1/(.std^2).
%       in the case of logarithmic parametrisation this corresponds to a
%       normal prior with mean .mu and variance (.std)^2.
%   .exp ... cell array containing specific information about individual
%       experiments. must have the following fields
%       .PA_post_processing ... (optional, only for PA data) function
%       handle for post-processing of PA for e.g. normalization
%       .sigma_noise ... function of the mixed effect parameter yielding
%           the standard deviation in measurements
%       .sigma_time ... function of the mixed effect parameter yielding
%           the standard deviation in event time-points
%       .noise_model ... string for the noise model for measurements.
%           either 'normal' or 'lognormal'
%       .parameter_model ... distribution assumption for the random effect.
%           either 'normal' or 'lognormal'
%       .model ... function handle for simulation of the model with arguments
%               t ... time
%               phi ... mixed effect parameter
%               kappa ... experimental condtion
%               option_model ... struct which can carry additional options such
%                   as the number of required sensitivities
%           the function handle should return the following object
%               sol ... solution struct with the following fields
%                   (depending on the value of option_model.sensi)
%                   for option_model.sensi >= 0
%                   .status ... >= 0 for successful simulation
%                   .y ... model output for observable, matrix [time,measurement]
%                   .root ... (optional) model output for event, matrix [time,event]
%                   .rootval ... (optional)  model output for rootfunction, matrix [time,event]
%                   for option_model.sensi >= 1
%                   .sy ... sensitivity for observable, matrix [time,measurement]
%                   .sroot ... (optional) sensitivity for event, matrix [time,event]
%                   .srootval ... (optional) sensitivity for rootfunction, matrix [time,event]
%                   for option_model.sensi >= 2 (optional)
%                   .s2y ... second order sensitivity for observable, matrix [time,measurement]
%                   .s2root ... (optional) second order sensitivity for event, matrix [time,event]
%                   .s2rootval ... (optional) second order sensitivity for rootfunction, matrix [time,event]
%      .fh ...figure handle for plotting of simulation
%      .fp ...figure handle for plotting of random effect distribution
%      .fl ...figure handle for plotting of likelihood contributions
%      .plot ... function handle for plotting of simulation resutls with arguments
%          Data ... data
%          Sim ... simulation
%          fh ... figure handle to the figure in which the function
%               should plot
%      and the following fields which contain figure handles whith mixed
%      effect parameter as argument of the respective variable
%      .dsigma_noisedphi
%      .ddsigma_noisedphidphi
%      .dddsigma_noisedphidphidphi
%      .ddddsigma_noisedphidphidphidphi (only for .integration==1)
%      .dsigma_timedphi
%      .ddsigma_timedphidphi
%      .dddsigma_timedphidphidphi
%      .ddddsigma_timedphidphidphidphi (only for .integration==1)
%      and the following fields which contain figure handles whith optimization
%      parameter as argument of the respective variable
%      .beta ... common effect parameter
%      .delta ... parametrization of covariance matrix for random effect
%      .dbetadxi
%      .ddeltadxi
%      .ddbetadxidxdi
%      .dddeltadxidxdi
%      and the following fields which contain figure handles whith random
%      and common effect parameter as argument
%      .phi ... mixed effect parameter @(b,beta)
%      .dphidbeta
%      .dphidb
%      .ddphidbetadbeta
%      .ddphidbdbeta
%      .ddphidbdb
%  options ... option struct with the following options
%      .tau_update ... minimum number of second which must pass before the
%      plots are updated
%      .plot ... flag whether the function should plot either
%          0 ... no plots
%          1 ... all plots (default)
%  extract_flag ... flag indicating whether the values of random effect
%      parameters are to be extracted (only for SCTL data)
%          0 ... no extraction (default)
%          1 ... extraction
%
%
% Outputs:
% ========
%      P ... cell array with field
%          .SCTL if the corresponding Data cell had a .SCTL field. this has
%          field has subfields
%              .bhat random effect parameter
%              .dbdxi gradient of random effect parameter wrt optimization
%              parameter
%              .ddbdxidxi hessian of random effect parameter wrt optimization
%              parameter
%      (otherwise)
%      SP.B ... location of sigma-points
%
% 2015/04/14 Fabian Froehlich

function varargout = validate_logL_CE_w_grad_2(varargin)

%% Load old values
persistent tau
persistent P_old
persistent logL_old

if isempty(tau)
    tau = clock;
end
if isempty(logL_old)
    logL_old = -inf;
end

%% Initialization
xi = varargin{1};
Data = varargin{2};
Model = varargin{3};

% Options
options.tau_update = 0;
options.plot = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

if nargin >= 5
    extract_flag = varargin{5};
else
    extract_flag = false;
end

nderiv = nargout;

% Plot options
if (etime(clock,tau) > options.tau_update) && (options.plot == 1)
    options.plot = 30;
    tau = clock;
else
    options.plot = 0;
end

%% Evaluation of likelihood function
% Initialization
logL = 0;
if nderiv >= 2
    dlogLdxi = zeros(length(xi),1);
    if nderiv >= 3
        ddlogLdxidxi = zeros(length(xi));
    end
end

% definition of possible datatypes
data_type = {'SCTL','SCSH','SCTLstat','PA'};

% Loop: Experiments/Experimental Conditions
for s = 1:length(Data)
    
    %% Assignment of global variables
    type_D = Model.type_D;
    
    %% Construct fixed effects and covariance matrix
    beta = Model.exp{s}.beta(xi);
    delta = Model.exp{s}.delta(xi);
    
    [D,invD] = xi2D(delta,type_D);
    
    % debugging:
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,1,3)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,3,5)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,2,4)
    % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,type_D),1e-4,4,6)
    
    
    %% Construction of time vector
    t_s = [];
    for dtype = 1:length(data_type)
        if isfield(Data{s},data_type{dtype})
            t_s = union(eval(['Data{s}.' data_type{dtype} '.time']),t_s);
        end
    end
    
    %% Single cell time-lapse data - Individuals
    if isfield(Data{s},'SCTL')
        % Reset values of Data likelihood and parameter likelihood
        
        % Evaluation of time index set
        [~,ind_time] = ismember(Data{s}.SCTL.time,t_s);
        
        % Initialization
        % measurements
        Sim_SCTL.Y = nan(size(Data{s}.SCTL.Y));
        % events
        if(~isfield(Data{s}.SCTL,'T'))
            Data{s}.SCTL.T = zeros(0,1,size(Data{s}.SCTL.Y,3));
        end
        Sim_SCTL.T = nan(size(Data{s}.SCTL.T));
        Sim_SCTL.R = nan(size(Data{s}.SCTL.T));
        
        % set default scaling
        if(~isfield(Model,'SCTLscale'))
            Model.SCTLscale = 1;
        end
        
        % load values from previous evaluation as initialisation
        if logL_old == -inf
            bhat_0 = zeros(length(Model.exp{s}.ind_b),size(Data{s}.SCTL.Y,3));
        else
            bhat_0 = P_old{s}.SCTL.bhat;
        end
        
        % Loop: Indiviudal cells
        
        dbetadxi = Model.exp{s}.dbetadxi(xi);
        ddeltadxi = Model.exp{s}.ddeltadxi(xi);
        ddbetadxidxi = Model.exp{s}.ddbetadxidxi(xi);
        dddeltadxidxi = Model.exp{s}.dddeltadxidxi(xi);
        
        logLi_D = zeros(1,size(Data{s}.SCTL.Y,3));
        logLi_T = zeros(1,size(Data{s}.SCTL.Y,3));
        logLi_b = zeros(1,size(Data{s}.SCTL.Y,3));
        logLi_I = zeros(1,size(Data{s}.SCTL.Y,3));
        
        for i = 1:size(Data{s}.SCTL.Y,3)
            
            % Load single-cell data
            Ym_si = Data{s}.SCTL.Y(ind_time,:,i);
            ind_y = find(~isnan(Ym_si));
            
            Tm_si = Data{s}.SCTL.T(:,:,i);
            ind_t = find(~isnan(Tm_si));
            
            bhat_si0 = bhat_0(:,i);
            
            %% Estimation of single cell random effects
            % Higher order derivatives of the objective function for single cell parameters
            % here G is the Hessian of the objective function and bhat_si is the optimum of the objective function
            % with respect to b
            % F_diff and b_diff determine how many derivatives of the objective function and of the optimum need to
            % be computed
            F_diff = 0;
            b_diff = 0;
            [bhat_si,G] ...
                = optimize_SCTL_si(Model,Data,bhat_si0,beta,delta,type_D,t_s,Ym_si,Tm_si,ind_y,ind_t,F_diff,b_diff,s);
            
            
            % Store bhat
            bhat(:,i) = bhat_si;
            
            % Construct single-cell parameter
            phi_si = Model.exp{s}.phi(beta,bhat_si);
            
            % Simulate model and compute derivatives
            [Y_si,T_si,R_si] = simulate_trajectory(t_s,phi_si,Model,Data{s}.condition,s,ind_t,ind_y);

            % Construct sigma
            [Sigma_noise_si] = build_sigma_noise(phi_si,Ym_si,s,Model,ind_y);
            [Sigma_time_si] = build_sigma_time(phi_si,Tm_si,s,Model,ind_t);

            
            %% Evaluation of likelihood and likelihood gradient
            
            % this is part accounts for the noise model
            % J_D = log(p(Y(b,beta)|D))
            switch(Model.exp{s}.noise_model)
                case 'normal'
                    J_D = normal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
                case 'lognormal'
                    J_D = lognormal_noise(Y_si,Ym_si,Sigma_noise_si,ind_y);
            end
            
            % this is part accounts for the event model
            % J_D = log(p(Y(b,beta)|D))
            switch(Model.exp{s}.noise_model)
                case 'normal'
                    J_T = normal_time(T_si,Tm_si,R_si,Sigma_time_si,ind_t);
            end
            
            % this part accounts for the parameter model
            % J_b = log(p(b_si|delta))
            switch(Model.exp{s}.parameter_model)
                case 'normal'
                    J_b = normal_param(bhat_si,delta,type_D);
                case 'lognormal'
                    J_b = lognormal_param(bhat_si,delta,type_D);
            end
            
            logLi_D(1,i) =  - J_D;
            logLi_T(1,i)  = - J_T;
            logLi_b(1,i)  = - J_b;
            
            if(Model.integration)
                % laplace approximation
                logLi_I(1,i) = - 0.5*log(det(G));
            end
            
            Y_si_tmp{i} = Y_si;
            T_si_tmp{i} = T_si;
            R_si_tmp{i} = R_si;
            
        end
        
        for k = 1:size(Data{s}.SCTL.Y,3)
            
            Ym_si = Data{s}.SCTL.Y(ind_time,:,k);
            idx_y = find(~isnan(Ym_si));
            
            Tm_si = Data{s}.SCTL.T(:,:,k);
            idx_t = find(~isnan(Tm_si));
            
            % Assignment of simulation results
            [I,J] = ind2sub([size(Sim_SCTL.Y,1),size(Sim_SCTL.Y,2)],idx_y);
            for i_ind = 1:length(I)
                Y_si = Y_si_tmp{k};
                Sim_SCTL.Y(I(i_ind),J(i_ind),k) = Y_si(i_ind);
            end
            [I,J] = ind2sub([size(Sim_SCTL.Y,1),size(Sim_SCTL.Y,2)],idx_t);
            for i_ind = 1:length(I)
                T_si = T_si_tmp{k};
                R_si = R_si_tmp{k};
                Sim_SCTL.T(I(i_ind),J(i_ind),k) = T_si(i_ind);
                Sim_SCTL.R(I(i_ind),J(i_ind),k) = R_si(i_ind);
            end
        end
        
        P{s}.SCTL.bhat = bhat;
        
        logL = logL + Model.SCTLscale*sum(logLi_D + logLi_T + logLi_b + logLi_I,2);
        
        if(~isfield(Model,'shr_fun'))
            if Model.penalty > 3
                error('Shrinkage Function is missing! please provide Model.shr_fun')
            else
                Model.shr_fun = 0;
            end
        end
        
        if(Model.penalty)
            % parameter penalization terms
            % logL_s = log(p(mu_S,S_s|b_s,D))
            logL_s = penal_param(P{s}.SCTL.bhat,delta,type_D,Model.penalty,Model.shr_fun);

            logL = logL + logL_s;
        end
        
        % Visulization
        if options.plot
            
            % Visualisation of single cell parameters
            figure(Model.exp{s}.fp)
            clf
            b_s = P{s}.SCTL.bhat;
            n_b = size(b_s,1);
            
            for j = 1:n_b
                subplot(ceil((n_b+1)/4),4,j+1)
                xx = linspace(-5*sqrt(D(j,j)),5*sqrt(D(j,j)),100);
                %nhist(P{s}.SCTL.bhat(j,:),'pdf','noerror');
                hold on
                plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','LineWidth',2)
                ecdf = zeros(length(xx),1);
                for k = 1:length(xx)
                    ecdf(k) = sum(b_s(j,:)<xx(k))/length(b_s(j,:));
                end
                plot(xx,ecdf,'--r','LineWidth',2)
                
                
                if(j==1)
                    
                end
                xlim([-5*sqrt(D(j,j)),5*sqrt(D(j,j))])
                ylim([0,1.1])
                %xlabel(char(Model.sym.b(Model.exp{s}.ind_b(j))));
                ylabel('cdf')
                box on
            end
            subplot(ceil(n_b+1/4),4,1,'Visible','off')
            hold on
            plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','Visible','off')
            plot(xx,ecdf,'--r','LineWidth',2,'Visible','off')
            
            
            legend('cdf of single cell Parameters','cdf of populaton Parameters')
            
            % Visualisation of likelihood contribution
            
            figure(Model.exp{s}.fl)
            if(Model.penalty)
                if(Model.integration)
                    bar([logLi_D;logLi_T;logLi_b;logLi_I;repmat(logL_s,[1,size(Data{s}.SCTL.Y,3)])],'stacked')
                    set(gca,'XTickLabel',{'Data','Event','Par','Int','Pen'})
                else
                    bar([logLi_D;logLi_T;logLi_b;repmat(logL_s,[1,size(Data{s}.SCTL.Y,3)])],'stacked')
                    set(gca,'XTickLabel',{'Data','Event','Par','Pen'})
                end
            else
                if(Model.integration)
                    bar([logLi_D;logLi_T;logLi_b;logLi_I],'stacked')
                    set(gca,'XTickLabel',{'Data','Event','Par','Int'})
                else
                    bar([logLi_D;logLi_T;logLi_b],'stacked')
                    set(gca,'XTickLabel',{'Data','Event','Par'})
                end
            end
            ylabel('log-likelihood')
            title('likelihood contribution')
            
            
            
            % Visualisation of data and fit
            Model.exp{s}.plot(Data{s},Sim_SCTL,Model.exp{s}.fh);
        end
    end
    
    %% Single cell time-lapse data - Statistics
    if isfield(Data{s},'SCTLstat')
        % Simulation using sigma points
        op_SP.nderiv = nderiv;
        op_SP.req = [0,0,0,1,1,1,0];
        if(nderiv>0)
            dbetadxi = Model.exp{s}.dbetadxi(xi);
            ddeltadxi = Model.exp{s}.ddeltadxi(xi);
            SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition),beta,D,Model.exp{s},op_SP,dDddelta,ddeltadxi,dbetadxi);
        else
            SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition),beta,D,Model.exp{s},op_SP);
        end
        
        % Evaluation of likelihood, likelihood gradient and hessian
        
        % Mean
        logL_mz = - 0.5*sum(nansum(((Data{s}.SCTLstat.mz - SP.mz)./Data{s}.SCTLstat.Sigma_mz).^2,1),2);
        if nderiv >= 2
            dlogL_mzdxi = permute(nansum(bsxfun(@times,(Data{s}.SCTLstat.mz - SP.mz)./Data{s}.SCTLstat.Sigma_mz.^2,SP.dmzdxi),1),[2,1]);
            if nderiv >= 3
                wdmz_SP = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_mz,SP.dmzdxi);
                %                     wdmz_SP = reshape(wdmz_SP,[numel(SP.mz),size(SP.dmdxizdxi,3)]);
                ddlogL_mzdxi2 = -wdmz_SP'*wdmz_SP;
            end
        end
        
        
        % Covariance
        logL_Cz = - 0.5*sum(nansum(nansum(((Data{s}.SCTLstat.Cz - SP.Cz)./Data{s}.SCTLstat.Sigma_Cz).^2,1),2),3);
        if nderiv >= 2
            dlogL_Czdxi = squeeze(nansum(nansum(bsxfun(@times,(Data{s}.SCTLstat.Cz - SP.Cz)./Data{s}.SCTLstat.Sigma_Cz.^2,SP.dCzdxi),1),2));
            if nderiv >= 3
                wdCzdxi = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_Cz,SP.dCzdxi);
                wdCzdxi = reshape(wdCzdxi,[numel(SP.Cz),size(SP.dCzdxi,3)]);
                ddlogL_Czdxi2 = -wdCzdxi'*wdCzdxi;
            end
        end
        
        % Summation
        logL = logL + logL_mz + logL_Cz;
        if nderiv >=2
            dlogLdxi = dlogLdxi + dlogL_mzdxi + dlogL_Czdxi;
            if nderiv >=3
                ddlogLdxidxi = ddlogLdxidxi + ddlogL_mzdxi2 + ddlogL_Czdxi2;
            end
        end
        
        % Visulization
        if options.plot
            Sim_SCTLstat.mz = SP.mz;
            Sim_SCTLstat.Cz = SP.Cz;
            Model.exp{s}.plot(Data{s},Sim_SCTLstat,Model.exp{s}.fh);
        end
        
    end
    
    %% Single cell snapshot data
    if isfield(Data{s},'SCSH')
        % Simulation using sigma points
        op_SP.nderiv = nderiv;
        op_SP.req = [1,1,0,0,0,1,0];
        op_SP.type_D = Model.type_D;
        SP = testSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCSH.time,phi,Data{s}.condition),xi,Model.exp{s},op_SP);
        
        % Evaluation of likelihood, likelihood gradient and hessian
        % Mean
        logL_m = - 0.5*nansum(nansum(((Data{s}.SCSH.m - SP.my)./Data{s}.SCSH.Sigma_m).^2,1),2);
        if nderiv >= 2
            dlogL_mdxi = squeeze(nansum(nansum(bsxfun(@times,(Data{s}.SCSH.m - SP.my)./Data{s}.SCSH.Sigma_m.^2,SP.dmydxi),1),2));
            if nderiv >= 3
                wdmdxi = bsxfun(@times,1./Data{s}.SCSH.Sigma_m,SP.dmydxi);
                wdmdxi = reshape(wdmdxi,[numel(SP.my),size(SP.dmydxi,3)]);
                ddlogL_mdxi2 = -wdmdxi'*wdmdxi;
            end
        end
        
        % Covariance
        logL_C = - 0.5*sum(nansum(nansum(((Data{s}.SCSH.C - SP.Cy)./Data{s}.SCSH.Sigma_C).^2,1),2),3);
        if nderiv >= 2
            dlogL_Cdxi = squeeze(nansum(nansum(nansum(bsxfun(@times,(Data{s}.SCSH.C - SP.Cy)./Data{s}.SCSH.Sigma_C.^2,SP.dCydxi),1),2),3));
            if nderiv >= 3
                wdCdxi = bsxfun(@times,1./Data{s}.SCSH.Sigma_C,SP.dCydxi);
                wdCdxi = reshape(wdCdxi,[numel(SP.Cy),size(SP.dCydxi,4)]);
                ddlogL_Cdxi2 = -wdCdxi'*wdCdxi;
            end
        end
        
        % Summation
        logL = logL + logL_m + logL_C;
        if nderiv >= 2
            dlogLdxi = dlogLdxi + dlogL_mdxi + dlogL_Cdxi;
            if nderiv >= 3
                ddlogLdxidxi = ddlogLdxidxi + ddlogL_mdxi2 + ddlogL_Cdxi2;
            end
        end
        
        % Visulization
        if options.plot
            Sim_SCSH.m = SP.my;
            Sim_SCSH.C = SP.Cy;
            Model.exp{s}.plot(Data{s},Sim_SCSH,Model.exp{s}.fh);
        end
        
    end
    
    %% Population average data
    if isfield(Data{s},'PA')
        
        % Simulation using sigma points
        op_SP.nderiv = nderiv;
        op_SP.req = [1,0,0,0,0,1,0];
        if(nderiv>0)
            dbetadxi = Model.exp{s}.dbetadxi(xi);
            ddeltadxi = Model.exp{s}.ddeltadxi(xi);
            SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.PA.time,phi,Data{s}.condition),beta,D,Model.exp{s},op_SP,dDddelta,ddeltadxi,dbetadxi);
        else
            SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.PA.time,phi,Data{s}.condition),beta,D,Model.exp{s},op_SP);
        end
        
        % Post-processing of population average data
        if isfield(Model.exp{s},'PA_post_processing')
            if(nderiv==1)
                SP.dmydxi = zeros([size(SP.my) size(xi,1)]);
            end
            [SP.my,SP.dmydxi] = Model.exp{s}.PA_post_processing(SP.my,SP.dmydxi);
        end
        
        
        % Evaluation of likelihood, likelihood gradient and hessian
        logL_m = - 0.5*nansum(nansum(((Data{s}.PA.m - SP.my)./Data{s}.PA.Sigma_m).^2,1),2);
        if nderiv >= 2
            dlogL_mdxi = squeeze(nansum(nansum(bsxfun(@times,(Data{s}.PA.m - SP.my)./Data{s}.PA.Sigma_m.^2,SP.dmydxi),1),2));
            if nderiv >= 3
                wdmdxi = bsxfun(@times,1./Data{s}.PA.Sigma_m,SP.dmydxi);
                wdmdxi = reshape(wdmdxi,[numel(SP.my),size(SP.dmydxi,3)]);
                ddlogL_mdxi2 = -wdmdxi'*wdmdxi;
            end
        end
        
        % Summation
        logL = logL + logL_m;
        if nderiv >= 2
            dlogLdxi = dlogLdxi + dlogL_mdxi;
            if nderiv >= 3
                ddlogLdxidxi = ddlogLdxidxi + ddlogL_mdxi2;
            end
        end
        
        % Visulization
        if options.plot
            Sim_PA.m = SP.m;
            Model.exp{s}.plot(Data{s},Sim_PA,Model.exp{s}.fh);
        end
    end
    
end

%% Output

if extract_flag
    if isfield(Data{s},'SCTL')
        varargout{1} = P;
    elseif any([isfield(Data{s},'SCTLstat'),isfield(Data{s},'PA'),isfield(Data{s},'SCSH')])
        varargout{1} = SP.B;
    end
    return
end

%% Prior

if isfield(Model,'prior')
    if(iscell(Model.prior))
        if(length(Model.prior) <= length(xi))
            for ixi = 1:length(Model.prior)
                if(isfield(Model.prior{ixi},'mu') && isfield(Model.prior{ixi},'std'))
                    if nderiv >= 1
                        % One output
                        logL =  logL + 0.5*((xi(ixi)-Model.prior{ixi}.mu)/Model.prior{ixi}.std)^2;
                        if nderiv >= 2
                            % Two outputs
                            dlogLdxi(ixi) =  dlogLdxi(ixi) + ((xi(ixi)-Model.prior{ixi}.mu)/Model.prior{ixi}.std^2);
                            if nderiv >= 3
                                % Two outputs
                                ddlogLdxidxi(ixi,ixi) =  ddlogLdxidxi(ixi,ixi) + 1/Model.prior{ixi}.std^2;
                            end
                        end
                    end
                end
            end
        else
            error('length of Model.prior must agree with length of optimization parameter')
        end
    else
        error('Model.prior must be a cell array')
    end
end


%%

if nderiv >= 1
    % One output
    varargout{1} =  logL;
    if nderiv >= 2
        % Two outputs
        varargout{2} =  dlogLdxi;
        if nderiv >= 3
            % Two outputs
            varargout{3} =  ddlogLdxidxi;
        end
    end
end

end