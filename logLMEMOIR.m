% logL_CE_w_grad_2 is the main likelihood function of the MEMToolbox and
% supports classical mixed effect models aswell as sigma-point based
% moment-approximation schemes to the solution of population balance
% equations. The likelihood function allows for the parallel estimation of
% parameters from multiple different experimental setups/conditions
%
% USAGE:
% ======
% [logL,dlogLdxi,ddlogLdxidxi] = logL_CE_w_grad_2(xi,Data,Model,options)
% P = logL_CE_w_grad_2(xi,Data,Model,options,1)
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
%       'matrix-logarithm' for full matrix with log. parametrisation
%   .integration ... flag indicating whether integration for classical
%       mixed effect models via laplace approximation should be applied.
%       only applicable for SCTL data.
%   .penalty ... flag indicating whether additional penalty terms for
%       synchronisation of parameters across experiments should be applied.
%       only applicable for SCTL data.
%   .prior ... cell array containing prior information for the optimization
%       parameters. the length of the cell array must not exceed the
%       length of the optimization parameter. the prior is applied when
%       there exist fields .mu and .std. This imposes a normal prior with
%       mean .mu and variance (.std)^2.
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
%      .plot ... function handle for plotting of simulation resutls with arguments
%          Data ... data
%          Sim ... simulation
%          s ... experiment index
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
%      .ms_iter ... number of function call with decreasing likelihood
%      value between multistarts for the inner optimisation problem for SCTL
%      data, for optimisation this is typically the number of iterations
%      between multistarts
%          0 ... multistart every iteration
%          X ... multistart every X iterations (default = 10)
%        Inf ... Only one multistart in the beginning
%  extract_flag ... flag indicating whether the values of random effect
%      parameters are to be extracted (only for SCTL data)
%          0 ... no extraction (default)
%          1 ... extraction
%
%
% Outputs:
% ========
%  extract_flag == 0
%      logL ... logLikelihood value
%      dlogLdxi ... gradient of logLikelihood
%      ddlogLdxidxdi ... hessian of logLikelihood
%  extract_flag == 1
%      (SCTL)
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

function varargout = logLMEMOIR(varargin)
   %% Load old values
   persistent tau
   persistent P_old
   persistent logL_old
   tic;
   if isempty(tau)
      tau = clock;
   end
   
   %% Initialization
   xi = varargin{1};
   Data = varargin{2};
   Model = varargin{3};
   
   % Options
   options.tau_update = 0;
   options.plot = 0;
   options.ms_iter = 10;
   options.events = 1;
   options.rescaleSCTL = 0;
   options.optimal_sigma = 1;
   if nargin >= 4
      if(isstruct(varargin{4}))
         options = setdefault(varargin{4},options);
      end
   end
   if nargin >= 5
      extract_flag = varargin{5};
   else
      extract_flag = false;
   end
   if nargin >= 6
      P_old = varargin{6};
      logL_old = -Inf;
   end
   options.nderiv = max(nargout-1,0);
   
   
   % initialise storage
   if(isempty(logL_old))
      logL_old = -Inf;
      for s = 1:length(Data)
         P_old{s}.n_store = 0;
         P_old{s}.xi = zeros(size(xi));
         if isfield(Data{s},'SCTL')
            P_old{s}.SCTL.bhat = zeros(length(Model.exp{s}.ind_b),size(Data{s}.SCTL.Y,3));
            P_old{s}.SCTL.dbdxi = zeros(length(Model.exp{s}.ind_b),length(xi),size(Data{s}.SCTL.Y,3));
         end
      end
   end
   P = P_old;
   
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
   if options.nderiv >= 1
      dlogLdxi = zeros(length(xi),1);
      if options.nderiv >= 2
         ddlogLdxidxi = zeros(length(xi));
      end
   end
   
   options.type_D = Model.type_D;
   options.integration = Model.integration;
   
   % Loop: Experiments/Experimental Conditions
   for s = 1:length(Data)
      
      %% Single cell time-lapse data - Individuals
      if isfield(Data{s},'SCTL')
         switch(options.nderiv)
            case 0
               [P,logL_sc] = logL_SCTL(xi, Model.exp{s}, Data{s}, s, options, P);
            case 1
               [P,logL_sc,dlogL_scdxi] = logL_SCTL(xi, Model.exp{s}, Data{s}, s, options, P);
               % [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logL_SCTL(xi, Model.exp{s}, Data{s}, s, options, P),1e-3,2,3)
            case 2
               [P,logL_sc,dlogL_scdxi,ddlogL_scdxi2] = logL_SCTL(xi, Model.exp{s}, Data{s}, s, options ,P);
         end
         
         if(options.rescaleSCTL)
            scaling = 1/size(Data{s}.SCTL.Y,3);
         else
            scaling = 1;
         end
         
         logL = logL + scaling*sum(logL_sc);
         if options.nderiv >= 1
            dlogLdxi = dlogLdxi + transpose(scaling*sum(dlogL_scdxi));
            if options.nderiv >= 2
               ddlogLdxidxi = ddlogLdxidxi + permute(scaling*sum(ddlogL_scdxi2),[2,3,1]);
            end
         end
      end
      
      %% Single cell time-lapse data - Statistics
      if isfield(Data{s},'SCTLstat')
         % Simulation using sigma points
         op_SP.nderiv = options.nderiv;
         op_SP.req = [0,0,0,1,1,1,0];
         op_SP.type_D = Model.type_D;
         % op_SP.approx = 'sample';
%          if(extract_flag)
%             SP = testSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition,'lin'),xi,Model.exp{s},op_SP);
%          else
            SP = getSigmaPointApp(@(phi) simulateForSP(Model.exp{s}.model,Data{s}.SCTLstat.time,phi,Data{s}.condition,'lin'),xi,Model.exp{s},op_SP);
%          end
         
         % Evaluation of likelihood, likelihood gradient and hessian
         
         % Mean
         logL_mz = - 0.5*sum(nansum(((Data{s}.SCTLstat.mz - SP.mz)./Data{s}.SCTLstat.Sigma_mz).^2,1),2);
         if options.nderiv >= 1
            dlogL_mzdxi = permute(nansum(bsxfun(@times,(Data{s}.SCTLstat.mz - SP.mz)./Data{s}.SCTLstat.Sigma_mz.^2,SP.dmzdxi),1),[2,1]);
            if options.nderiv >= 2
               wdmz_SP = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_mz,SP.dmzdxi);
               %                     wdmz_SP = reshape(wdmz_SP,[numel(SP.mz),size(SP.dmdxizdxi,3)]);
               ddlogL_mzdxi2 = -wdmz_SP'*wdmz_SP;
            end
         end
         
         
         % Covariance
         logL_Cz = - 0.5*sum(nansum(nansum(((Data{s}.SCTLstat.Cz - SP.Cz)./Data{s}.SCTLstat.Sigma_Cz).^2,1),2),3);
         if options.nderiv >= 1
            dlogL_Czdxi = squeeze(nansum(nansum(bsxfun(@times,(Data{s}.SCTLstat.Cz - SP.Cz)./Data{s}.SCTLstat.Sigma_Cz.^2,SP.dCzdxi),1),2));
            if options.nderiv >= 2
               wdCzdxi = bsxfun(@times,1./Data{s}.SCTLstat.Sigma_Cz,SP.dCzdxi);
               wdCzdxi = reshape(wdCzdxi,[numel(SP.Cz),size(SP.dCzdxi,3)]);
               ddlogL_Czdxi2 = -wdCzdxi'*wdCzdxi;
            end
         end
         
         % Summation
         logL = logL + logL_mz + logL_Cz;
%          disp([num2str(s) '     ' num2str(logL) '     ' num2str(logL_mz) '     ' num2str(logL_Cz)])
         if options.nderiv >= 1
            dlogLdxi = dlogLdxi + dlogL_mzdxi + dlogL_Czdxi;
            if options.nderiv >= 2
               ddlogLdxidxi = ddlogLdxidxi + ddlogL_mzdxi2 + ddlogL_Czdxi2;
            end
         end
         
         % Visulization
         if options.plot
            Sim_SCTLstat.mz = SP.mz;
            Sim_SCTLstat.Cz = SP.Cz;
            Model.exp{s}.plot(Data{s},Sim_SCTLstat,s);
         end
         
         P{s}.SCTLstat.SP = SP;
         
      end
      
      %% Single cell snapshot data
      if isfield(Data{s},'SCSH')
         
         switch(options.nderiv)
            case 0
               [SP,logL_m,logL_C] = logL_SCSH(xi, Model, Data, s, options);
            case 1
               [SP,logL_m,logL_C,dlogL_mdxi,dlogL_Cdxi] = logL_SCSH(xi, Model, Data, s, options);
            case 2
               [SP,logL_m,logL_C,dlogL_mdxi,dlogL_Cdxi,ddlogL_mdxi2,ddlogL_Cdxi2] = logL_SCSH(xi, Model, Data, s, options);
         end
         
         % Summation
         logL = logL + logL_m + logL_C;
         if options.nderiv >= 1
            dlogLdxi = dlogLdxi + dlogL_mdxi + dlogL_Cdxi;
            if options.nderiv >= 2
               ddlogLdxidxi = ddlogLdxidxi + ddlogL_mdxi2 + ddlogL_Cdxi2;
            end
         end
         
         P{s}.SCSH.SP = SP;
         
      end
      
      %% Population average data
      if isfield(Data{s},'PA')
         
         switch(options.nderiv)
            case 0
               [SP,logL_m] = logL_PA(xi, Model, Data, s, options);
            case 1
               [SP,logL_m,dlogL_mdxi] = logL_PA(xi, Model, Data, s, options);
               % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logL_PA(xi, Model, Data, s, options),1e-3,2,3)
            case 2
               [SP,logL_m,dlogL_mdxi,ddlogL_mdxi2] = logL_PA(xi, Model, Data, s, options);
         end
         
         % Summation
         logL = logL + logL_m;
         if options.nderiv >= 1
            dlogLdxi = dlogLdxi + dlogL_mdxi;
            if options.nderiv >= 2
               ddlogLdxidxi = ddlogLdxidxi + ddlogL_mdxi2;
            end
         end
         
         P{s}.PA.SP = SP;
      end
   end
   
   % updated stored value
   if(logL > logL_old)
      logL_old = logL;
      P_old = P;
      for s = 1:length(Data)
         P_old{s}.xi = xi;
         P_old{s}.n_store = P_old{s}.n_store + 1;
      end
   end
   
   
   %% Prior
   
   if isfield(Model,'prior')
      if(iscell(Model.prior))
         if(length(Model.prior) <= length(xi))
            for ixi = 1:length(Model.prior)
               if(isfield(Model.prior{ixi},'mu') && isfield(Model.prior{ixi},'std'))
                  if options.nderiv >= 1
                     % One output
                     logL =  logL - 0.5*((xi(ixi)-Model.prior{ixi}.mu)/Model.prior{ixi}.std)^2;
                     if options.nderiv >= 2
                        % Two outputs
                        dlogLdxi(ixi) =  dlogLdxi(ixi) - ((xi(ixi)-Model.prior{ixi}.mu)/Model.prior{ixi}.std^2);
                        if options.nderiv >= 3
                           % Two outputs
                           ddlogLdxidxi(ixi,ixi) =  dsdxi(ixi,ixi) - 1/Model.prior{ixi}.std^2;
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
   
   
   %% OUTPUT
   
   if extract_flag
      varargout{1} = P;
      return
   end
   
   if nargout >= 1
      % One output
      varargout{1} =  logL;
      if nargout >= 2
         % Two outputs
         varargout{2} =  dlogLdxi;
         if nargout >= 3
            % Two outputs
            varargout{3} =  ddlogLdxidxi;
         end
      end
   end
end