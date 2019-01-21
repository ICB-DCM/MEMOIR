% complete_model uses the symbolic definition of the  parametrisation in
% the input model struct to generate m-files which can be used to evaluate
% said symbolic expressions
%
% USAGE:
% ======
% MODEL = complete_model(MODEL,S)
%
% INPUTS:
% =======
% Model ... model struct encapsulating the model definition for a MEM
%   .sym ... contains symbolic definition of the overall model
%       .xi ... are the parameter wich are optimised, this usually consist
%       of common effecmts, the parametrisation of the random effects
%       covariance matrix and the parametrisation of the noise parameters
%       .phi ... is are the mixed effect parametrisation as function of
%       common effects beta and random effects b
%       .beta ...  is the parametrisation of common effects as function of
%       xi
%       .b ... is the parametrisation of random effects
%       .delta ... is the parametrisation of the covariance matrix. this
%       definition should be chosen in accordance to the definition of the
%       respective parametrisation given in Model.type_D
% S ... vector containing indexes of considered experiments
%
% Outputs:
% ========
% Model ... model struct encapsulating the model definition for a MEM
% problem
%   .exp{s} ... contains all information with respect to experiment number
%   s
%       .N ... number of single cells measured in the experiment
%       .sigma_noise ... single-cell noise level in this experiment (used 
%       in data generation and for fitting)
%       .sigma_mean ... noise level for mean measurements, used for fitting
%       .sigma_cov ... noise level of covariance and cross-covariance, used
%       .sigma_time ... noise level for event data (used in data
%       generation)
%       in fitting
%       .sigma_on ... flag indicating whether noise should be added during
%       data generation
%       .t ... vector of timepoints at which the system is observed
%       .ind_phi ... vector of indices of parameters which are active for
%       this experiment
%       .sym ... contains symbolic expression for the reduced parameters
%       for the respective experiment. moreover this struct will contain
%       the links to m-files for the evaluation of respective symbolic
%       expressions
%       .noise_model ... indicates the employed noise model for the
%       experiment
%       .parameter_model ... indicates the employed parameter model for
%       random effects
%       .fh ... figure handle for figure in which simulation for current
%       parameter values is compared against data
%       .fp ... figure handle for figure in which single cell parameters,
%       their empiric density aswell as their estimated density is plotted
%       .fl ... figure handle for figure in which the contribution of
%       individual terms to the objective function value is plotted
%       .plot ... function handle to function which generates the plot in
%       the figure for figure handle fh
% 2015/04/14 Fabian Froehlich


function Model = complete_model(Model,S)

% initialise model-loading flag
loadold = false;

% concatenate model name
switch(Model.type_D)
    case {'matrix-logarithm', 'givens-parametrization', 'Lie-generators'}
        filename = [Model.name '_full'];
    case 'diag-matrix-logarithm'
        filename = [Model.name '_diag'];
end
% generate path
[mdir,~,~]=fileparts(which(mfilename('full')));

%check for existence of directory
if(~exist(fullfile(mdir,'models',filename),'dir'))
    mkdir(fullfile(mdir,'models',filename));
end

% remove old paths and add new paths
addpath(fullfile(mdir,'models',filename));

try
    % check wheter the saved symbolic definition agrees with the current
    % one
    load(fullfile(mdir,'models',filename,'syms.mat'))
    f_xi = isequaln(Model.sym.xi,syms.xi);
    f_phi = isequaln(Model.sym.phi,syms.phi);
    f_beta = isequaln(Model.sym.beta,syms.beta);
    f_b = isequaln(Model.sym.b,syms.b);
    f_delta = isequaln(Model.sym.delta,syms.delta);
    for s = 1:length(Model.exp)
        
        if(isfield(Model.exp{s},'onlyPA'))
            onlyPA = Model.exp{s}.onlyPA;
        else
            onlyPA = false;
        end
        
        f_phiexp(s) = isequaln(Model.exp{s}.sym.phi,expsyms{s}.phi);
        f_sigma_noiseexp(s) = isequaln(Model.exp{s}.sym.sigma_noise,expsyms{s}.sigma_noise);
        f_sigma_timeexp(s) = isequaln(Model.exp{s}.sym.sigma_time,expsyms{s}.sigma_time);
        
        if(onlyPA)
            f_files(s) = all([exist(fullfile(mdir,'models',filename,['MEMbeta_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMdelta_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMdbetadxi_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMddeltadxi_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMsigma_noise_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMds_ndp_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMsigma_time_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMds_tdp_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMphi_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMdphidbeta_' filename '_' num2str(S(s))]),'file'),...
                exist(fullfile(mdir,'models',filename,['MEMdphidb_' filename '_' num2str(S(s))]),'file')]);
        else
        f_files(s) = all([exist(fullfile(mdir,'models',filename,['MEMbeta_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdelta_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdbetadxi_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddeltadxi_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddbetadxidxi_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdddeltadxidxi_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMsigma_noise_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMds_ndp_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdds_ndpdp_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddds_ndpdpdp_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMsigma_time_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMds_tdp_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdds_tdpdp_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddds_tdpdpdp_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMphi_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdphidbeta_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMdphidb_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddphidbdb_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddphidbetadbeta_' filename '_' num2str(S(s))]),'file'),...
            exist(fullfile(mdir,'models',filename,['MEMddphidbdbeta_' filename '_' num2str(S(s))]),'file')]);
        end
    end
    
    if(all([f_xi,f_phi,f_beta,f_b,f_delta,f_phiexp,f_sigma_noiseexp,f_sigma_timeexp,f_files]))
        %             if(Model.integration)
        %             eval(['Model.exp{s}.ddddsigma_timedphidphidphidphi = @MEMdddds_tdpdpdpdp_' filename '_' num2str(S(s)) ';']);
        %             end
        %             if(Model.integration)
        %             eval(['Model.exp{s}.ddddsigma_noisedphidphidphidphi = @MEMdddds_ndpdpdpdp_' filename '_' num2str(S(s)) ';']);
        %             end
        loadold = true;
        disp('Loading previous model definition files!');
        disp(['To regenerate model, abort and delete ' mdir 'models/' filename ]);
    end
catch
    
end

if(~loadold)
    % if we cannot load the old definition, we have to generate a new one
    disp(['Generating new model definition files!'])
    
    % save the symbolic definition as future reference
    syms = Model.sym;
    for s = 1:length(Model.exp)
        expsyms{s} = Model.exp{s}.sym;
    end
    save(fullfile(mdir,'models',filename,'syms.mat'),'syms','expsyms');
    
    % compute number of elements of xi and b
    n_xi = length(Model.sym.xi);
    n_b = length(Model.sym.b);
    
    % construct variance matrix parametrisation
    C = sym('C',[n_b,n_b]);
    switch(Model.type_D)
        case {'givens-parametrization','Lie-generators'}
            for j = 1:n_b
                C(j,j) = Model.sym.delta(j);
            end
            l = 1;
            for j = 1:n_b
                for k = 1:j-1
                    C(k,j) = Model.sym.delta(n_b + l);
                    C(j,k) = Model.sym.delta(n_b + l);
                    l = l + 1;
                end
            end
        case {'matrix-logarithm', 'cholesky-parametrization'}
            l = 1;
            for j = 1:n_b
                for k = 1:j
                    C(j,k) = Model.sym.delta(l);
                    C(k,j) = Model.sym.delta(l);
                    l = l+1;
                end
            end
        case 'diag-matrix-logarithm'
            C = diag(Model.sym.delta);
    end
    
    % loop over all experiments
    for s = 1:length(S)
        
        if(isfield(Model.exp{s},'onlyPA'))
            onlyPA = Model.exp{s}.onlyPA;
        else
            onlyPA = false;
        end
        
        % construct indices for reduced parameters
        Model.exp{s}.ind_beta = find(ismember(Model.sym.beta,symvar(Model.exp{s}.sym.phi)));
        if ~strcmp(Model.type_D, 'diag-matrix-logarithm')
            Model.exp{s}.ind_b = 1:length(Model.sym.b);
            Model.exp{s}.ind_delta = 1:length(Model.sym.delta);
        else
            Model.exp{s}.ind_b = find(ismember(Model.sym.b,symvar(Model.exp{s}.sym.phi)));
            Cs = C(Model.exp{s}.ind_b,Model.exp{s}.ind_b);
            Model.exp{s}.ind_delta =  find(ismember(Model.sym.delta,symvar(Cs)));
        end
        
        % constructe reduced parameters
        Model.exp{s}.sym.beta = Model.sym.beta(Model.exp{s}.ind_beta);
        Model.exp{s}.sym.b = Model.sym.b(Model.exp{s}.ind_b);
        Model.exp{s}.sym.delta = Model.sym.delta(Model.exp{s}.ind_delta);
        phi = Model.sym.phi(Model.exp{s}.ind_phi);
        
        % compute parameter length
        n_beta = length(Model.exp{s}.sym.beta);
        n_delta = length(Model.exp{s}.sym.delta);
        n_phi = length(Model.exp{s}.sym.phi);
        n_b = length(Model.exp{s}.sym.b);
        
        % generate m-files for parametrisation and respective derivatives
        % mfun is derived from 'matlabFunction' and largely follows the
        % same syntax but has some reduced functionality. we use mfun over
        % matlabFunction as matlabFunction does not adequately support
        % sparsity of symbolic variables which leads to a high
        % computational complexity even for relatively small models.
        
        % beta(xi) delta(xi)
        mfun(Model.exp{s}.sym.beta,'file',fullfile(mdir,'models',filename,['MEMbeta_' filename '_' num2str(S(s))]),'vars',{Model.sym.xi});
        eval(['Model.exp{s}.beta = @MEMbeta_' filename '_' num2str(S(s)) ';']);
        mfun(Model.exp{s}.sym.delta,'file',fullfile(mdir,'models',filename,['MEMdelta_' filename '_' num2str(S(s))]),'vars',{Model.sym.xi});
        eval(['Model.exp{s}.delta = @MEMdelta_' filename '_' num2str(S(s)) ';']);
        
        % dbetadxi
        Model.exp{s}.sym.dbetadxi = simplify(jacobian(Model.exp{s}.sym.beta,Model.sym.xi));
        mfun(Model.exp{s}.sym.dbetadxi,'file',fullfile(mdir,'models',filename,['MEMdbetadxi_' filename '_' num2str(S(s))]),'vars',{Model.sym.xi});
        eval(['Model.exp{s}.dbetadxi = @MEMdbetadxi_' filename '_' num2str(S(s)) ';']);
        
        % ddeltadxi
        Model.exp{s}.sym.ddeltadxi = simplify(jacobian(Model.exp{s}.sym.delta,Model.sym.xi));
        mfun(Model.exp{s}.sym.ddeltadxi,'file',fullfile(mdir,'models',filename,['MEMddeltadxi_' filename '_' num2str(S(s))]),'vars',{Model.sym.xi});
        eval(['Model.exp{s}.ddeltadxi = @MEMddeltadxi_' filename '_' num2str(S(s)) ';']);
        
        if(~onlyPA)
            % ddbetadxidxi
            Model.exp{s}.sym.ddbetadxidxi = sym(zeros(n_beta,n_xi,n_xi));
            for j = 1:n_beta
                Model.exp{s}.sym.ddbetadxidxi(j,:,:) = simplify(hessian(Model.exp{s}.sym.beta(j),Model.sym.xi));
            end
            mfun(Model.exp{s}.sym.ddbetadxidxi,'file',fullfile(mdir,'models',filename,['MEMddbetadxidxi_' filename '_' num2str(S(s))]),'vars',{Model.sym.xi});
            eval(['Model.exp{s}.ddbetadxidxi = @MEMddbetadxidxi_' filename '_' num2str(S(s)) ';']);
            
            % ddeltadxidxi
            Model.exp{s}.sym.dddeltadxidxi = sym(zeros(n_delta,n_xi,n_xi));
            for j = 1:n_delta
                Model.exp{s}.sym.dddeltadxidxi(j,:,:) = simplify(hessian(Model.exp{s}.sym.delta(j),Model.sym.xi));
            end
            mfun(Model.exp{s}.sym.dddeltadxidxi,'file',fullfile(mdir,'models',filename,['MEMdddeltadxidxi_' filename '_' num2str(S(s))]),'vars',{Model.sym.xi});
            eval(['Model.exp{s}.dddeltadxidxi = @MEMdddeltadxidxi_' filename '_' num2str(S(s)) ';']);
        end
        
        % sigma_noise(phi)
        mfun(Model.exp{s}.sym.sigma_noise,'file',fullfile(mdir,'models',filename,['MEMsigma_noise_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.sigma_noise = @MEMsigma_noise_' filename '_' num2str(S(s)) ';']);
        
        % dsigma_noisedphi
        Model.exp{s}.sym.dsigma_noisedphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
            for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
                Model.exp{s}.sym.dsigma_noisedphi(j,k,:) = jacobian(Model.exp{s}.sym.sigma_noise(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.dsigma_noisedphi,'file',fullfile(mdir,'models',filename,['MEMds_ndp_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.dsigma_noisedphi = @MEMds_ndp_' filename '_' num2str(S(s)) ';']);
        
        if(~onlyPA)
            % ddsigma_noisedphidphi
            Model.exp{s}.sym.ddsigma_noisedphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi,n_phi));
            for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
                for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
                    Model.exp{s}.sym.ddsigma_noisedphidphi(j,k,:,:) = hessian(Model.exp{s}.sym.sigma_noise(j,k),phi);
                end
            end
            mfun(Model.exp{s}.sym.ddsigma_noisedphidphi,'file',fullfile(mdir,'models',filename,['MEMdds_ndpdp_' filename '_' num2str(S(s))]),'vars',{phi});
            eval(['Model.exp{s}.ddsigma_noisedphidphi = @MEMdds_ndpdp_' filename '_' num2str(S(s)) ';']);
            
            % dddsigma_noisedphidphidphi
            Model.exp{s}.sym.dddsigma_noisedphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi,n_phi,n_phi));
            for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
                for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
                    for m = 1:n_phi
                        Model.exp{s}.sym.dddsigma_noisedphidphidphi(j,k,:,:,m) = diff(Model.exp{s}.sym.ddsigma_noisedphidphi(j,k,:,:),phi(m));
                    end
                end
            end
            mfun(Model.exp{s}.sym.dddsigma_noisedphidphidphi,'file',fullfile(mdir,'models',filename,['MEMddds_ndpdpdp_' filename '_' num2str(S(s))]),'vars',{phi});
            eval(['Model.exp{s}.dddsigma_noisedphidphidphi = @MEMddds_ndpdpdp_' filename '_' num2str(S(s)) ';']);
        end
        
        %         if(Model.integration)
        %             % ddddsigma_noisedphidphidphidphi --- currently disabled due to
        %             % high computational cost
        %             Model.exp{s}.sym.ddddsigma_noisedphidphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi,n_phi,n_phi,n_phi));
        % %             for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
        % %                 for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
        % %                     for m = 1:n_phi
        % %                         Model.exp{s}.sym.ddddsigma_noisedphidphidphidphi(j,k,:,:,:,m) = diff(Model.exp{s}.sym.dddsigma_noisedphidphidphi(j,k,:,:,:),phi(m));
        % %                     end
        % %                 end
        % %             end
        %             mfun(Model.exp{s}.sym.ddddsigma_noisedphidphidphidphi,'file',fullfile(mdir,'models',filename,['MEMdddds_ndpdpdpdp_' filename '_' num2str(S(s))]),'vars',{phi});
        %             eval(['Model.exp{s}.ddddsigma_noisedphidphidphidphi = @MEMdddds_ndpdpdpdp_' num2str(S(s)) ';']);
        %         end
        
        % sigma_mean(phi)
        mfun(Model.exp{s}.sym.sigma_mean,'file',fullfile(mdir,'models',filename,['MEMsigma_mean_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.sigma_mean = @MEMsigma_mean_' filename '_' num2str(S(s)) ';']);
        
        % dsigma_meandphi
        Model.exp{s}.sym.dsigma_meandphi = sym(zeros(size(Model.exp{s}.sym.sigma_mean,1),size(Model.exp{s}.sym.sigma_mean,2),n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_mean,1)
            for k = 1:size(Model.exp{s}.sym.sigma_mean,2)
                Model.exp{s}.sym.dsigma_meandphi(j,k,:) = jacobian(Model.exp{s}.sym.sigma_mean(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.dsigma_meandphi,'file',fullfile(mdir,'models',filename,['MEMds_mdp_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.dsigma_meandphi = @MEMds_mdp_' filename '_' num2str(S(s)) ';']);
        
        % sigma_cov(phi)
        mfun(Model.exp{s}.sym.sigma_cov,'file',fullfile(mdir,'models',filename,['MEMsigma_cov_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.sigma_cov = @MEMsigma_cov_' filename '_' num2str(S(s)) ';']);
        
        % dsigma_covdphi
        Model.exp{s}.sym.dsigma_covdphi = sym(zeros(size(Model.exp{s}.sym.sigma_cov,1),size(Model.exp{s}.sym.sigma_cov,2),n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_cov,1)
            for k = 1:size(Model.exp{s}.sym.sigma_cov,2)
                Model.exp{s}.sym.dsigma_covdphi(j,k,:) = jacobian(Model.exp{s}.sym.sigma_cov(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.dsigma_covdphi,'file',fullfile(mdir,'models',filename,['MEMds_cdp_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.dsigma_covdphi = @MEMds_cdp_' filename '_' num2str(S(s)) ';']);
        
        % sigma_time(phi)
        if(~isfield(Model.exp{s}.sym,'sigma_time'))
            Model.exp{s}.sym.sigma_time = sym.empty(0,1);
        end
        mfun(Model.exp{s}.sym.sigma_time,'file',fullfile(mdir,'models',filename,['MEMsigma_time_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.sigma_time = @MEMsigma_time_' filename '_' num2str(S(s)) ';']);
        
        % dsigma_timedphi
        Model.exp{s}.sym.dsigma_timedphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_time,1)
            for k = 1:size(Model.exp{s}.sym.sigma_time,2)
                Model.exp{s}.sym.dsigma_timedphi(j,k,:) = jacobian(Model.exp{s}.sym.sigma_time(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.dsigma_timedphi,'file',fullfile(mdir,'models',filename,['MEMds_tdp_' filename '_' num2str(S(s))]),'vars',{phi});
        eval(['Model.exp{s}.dsigma_timedphi = @MEMds_tdp_' filename '_' num2str(S(s)) ';']);
        
        if(~onlyPA)
            % ddsigma_timedphidphi
            Model.exp{s}.sym.ddsigma_timedphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi,n_phi));
            for j = 1:size(Model.exp{s}.sym.sigma_time,1)
                for k = 1:size(Model.exp{s}.sym.sigma_time,2)
                    Model.exp{s}.sym.ddsigma_timedphidphi(j,k,:,:) = hessian(Model.exp{s}.sym.sigma_time(j,k),phi);
                end
            end
            mfun(Model.exp{s}.sym.ddsigma_timedphidphi,'file',fullfile(mdir,'models',filename,['MEMdds_tdpdp_' filename '_' num2str(S(s))]),'vars',{phi});
            eval(['Model.exp{s}.ddsigma_timedphidphi = @MEMdds_tdpdp_' filename '_' num2str(S(s)) ';']);
            
            % dddsigma_timedphidphidphi
            Model.exp{s}.sym.dddsigma_timedphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi,n_phi,n_phi));
            for j = 1:size(Model.exp{s}.sym.sigma_time,1)
                for k = 1:size(Model.exp{s}.sym.sigma_time,2)
                    for m = 1:n_phi
                        Model.exp{s}.sym.dddsigma_timedphidphidphi(j,k,:,:,m) = diff(Model.exp{s}.sym.ddsigma_timedphidphi(j,k,:,:),phi(m));
                    end
                end
            end
            mfun(Model.exp{s}.sym.dddsigma_timedphidphidphi,'file',fullfile(mdir,'models',filename,['MEMddds_tdpdpdp_' filename '_' num2str(S(s))]),'vars',{phi});
            eval(['Model.exp{s}.dddsigma_timedphidphidphi = @MEMddds_tdpdpdp_' filename '_' num2str(S(s)) ';']);
        end
        %         if(Model.integration)
        %             % ddddsigma_timedphidphidphidphi --- currently disabled due to
        %             % high computational cost
        %             Model.exp{s}.sym.ddddsigma_timedphidphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi,n_phi,n_phi,n_phi));
        % %             for j = 1:size(Model.exp{s}.sym.sigma_time,1)
        % %                 for k = 1:size(Model.exp{s}.sym.sigma_time,2)
        % %                     for m = 1:n_phi
        % %                         Model.exp{s}.sym.ddddsigma_timedphidphidphidphi(j,k,:,:,:,m) = diff(Model.exp{s}.sym.dddsigma_timedphidphidphi(j,k,:,:,:),phi(m));
        % %                     end
        % %                 end
        % %             end
        %             mfun(Model.exp{s}.sym.ddddsigma_timedphidphidphidphi,'file',fullfile(mdir,'models',filename,['MEMdddds_tdpdpdpdp_' filename '_' num2str(S(s))]),'vars',{phi});
        %             eval(['Model.exp{s}.ddddsigma_timedphidphidphidphi = @MEMdddds_tdpdpdpdp_' filename '_' num2str(S(s)) ';']);
        %         end
        
        % phi
        mfun(Model.exp{s}.sym.phi,'file',fullfile(mdir,'models',filename,['MEMphi_' filename '_' num2str(S(s))]),'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.phi = @MEMphi_' filename '_' num2str(S(s)) ';']);
        
        % dphidbeta
        Model.exp{s}.sym.dphidbeta = simplify(jacobian(Model.exp{s}.sym.phi,Model.exp{s}.sym.beta));
        mfun(Model.exp{s}.sym.dphidbeta,'file',fullfile(mdir,'models',filename,['MEMdphidbeta_' filename '_' num2str(S(s))]),'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.dphidbeta = @MEMdphidbeta_' filename '_' num2str(S(s)) ';']);
        
        % dphidb
        Model.exp{s}.sym.dphidb = simplify(jacobian(Model.exp{s}.sym.phi,Model.exp{s}.sym.b));
        mfun(Model.exp{s}.sym.dphidb,'file',fullfile(mdir,'models',filename,['MEMdphidb_' filename '_' num2str(S(s))]),'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.dphidb = @MEMdphidb_' filename '_' num2str(S(s)) ';']);
        
        if(~onlyPA)
            % ddphidbdb
            Model.exp{s}.sym.ddphidbdb  = sym(zeros(n_phi,n_b,n_b));
            for j = 1:n_phi
                Model.exp{s}.sym.ddphidbdb(j,:,:) = simplify(hessian(Model.exp{s}.sym.phi(j),Model.exp{s}.sym.b));
            end
            mfun(Model.exp{s}.sym.ddphidbdb,'file',fullfile(mdir,'models',filename,['MEMddphidbdb_' filename '_' num2str(S(s))]),'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
            eval(['Model.exp{s}.ddphidbdb = @MEMddphidbdb_' filename '_' num2str(S(s)) ';']);
            
            % ddphidbetadbeta
            Model.exp{s}.sym.ddphidbetadbeta  = sym(zeros(n_phi,n_beta,n_beta));
            for j = 1:n_phi
                Model.exp{s}.sym.ddphidbetadbeta(j,:,:) = simplify(hessian(Model.exp{s}.sym.phi(j),Model.exp{s}.sym.beta));
            end
            mfun(Model.exp{s}.sym.ddphidbetadbeta,'file',fullfile(mdir,'models',filename,['MEMddphidbetadbeta_' filename '_' num2str(S(s))]),'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
            eval(['Model.exp{s}.ddphidbetadbeta = @MEMddphidbetadbeta_' filename '_' num2str(S(s)) ';']);
            
            % ddphidbdbeta
            Model.exp{s}.sym.ddphidbetadb  = sym(zeros(n_phi,n_beta,n_b));
            for j = 1:n_phi
                Model.exp{s}.sym.ddphidbdbeta(j,:,:) = simplify(jacobian(jacobian(Model.exp{s}.sym.phi(j),Model.exp{s}.sym.b),Model.exp{s}.sym.beta));
            end
            mfun(Model.exp{s}.sym.ddphidbdbeta,'file',fullfile(mdir,'models',filename,['MEMddphidbdbeta_' filename '_' num2str(S(s))]),'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
            eval(['Model.exp{s}.ddphidbdbeta = @MEMddphidbdbeta_' filename '_' num2str(S(s)) ';']);
        end
    end
    
else
    % if we can load the old definition, we just have to attach the m-files
    % to the model struct
    
    
    % save the symbolic definition as future reference
    syms = load(fullfile(mdir,'models',filename,'syms.mat'));
    Model.sym = syms.syms;
    
    % compute number of elements of xi and b
    n_xi = length(Model.sym.xi);
    n_b = length(Model.sym.b);
    
    
    % construct variance matrix parametrisation
    C = sym('C',[n_b,n_b]);
    switch(Model.type_D)
        case {'givens-parametrization','Lie-generators'}
            for j = 1:n_b
                C(j,j) = Model.sym.delta(j);
            end
            l = 1;
            for j = 1:n_b
                for k = 1:j-1
                    C(k,j) = Model.sym.delta(n_b + l);
                    C(j,k) = Model.sym.delta(n_b + l);
                    l = l + 1;
                end
            end
        case {'matrix-logarithm', 'cholesky-parametrization'}
            l = 1;
            for j = 1:n_b
                for k = 1:j
                    C(j,k) = Model.sym.delta(l);
                    C(k,j) = Model.sym.delta(l);
                    l = l+1;
                end
            end
        case 'diag-matrix-logarithm'
            C = diag(Model.sym.delta);
    end
    
    % loop over experiments
    for s = 1:length(Model.exp)
        
        if(isfield(Model.exp{s},'onlyPA'))
            onlyPA = Model.exp{s}.onlyPA;
        else
            onlyPA = false;
        end
        
        Model.exp{s}.ind_beta = find(ismember(Model.sym.beta,symvar(Model.exp{s}.sym.phi)));
        Model.exp{s}.ind_b = find(ismember(Model.sym.b,symvar(Model.exp{s}.sym.phi)));
        Cs = C(Model.exp{s}.ind_b,Model.exp{s}.ind_b);
        Model.exp{s}.ind_delta =  find(ismember(Model.sym.delta,symvar(Cs)));
        
        % construct reduced parameters
        Model.exp{s}.sym.beta = Model.sym.beta(Model.exp{s}.ind_beta);
        Model.exp{s}.sym.b = Model.sym.b(Model.exp{s}.ind_b);
        Model.exp{s}.sym.delta = Model.sym.delta(Model.exp{s}.ind_delta);
        
        eval(['Model.exp{s}.beta = @MEMbeta_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.delta = @MEMdelta_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dbetadxi = @MEMdbetadxi_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddeltadxi = @MEMddeltadxi_' filename '_' num2str(S(s)) ';']);
        if(~onlyPA)
            eval(['Model.exp{s}.ddbetadxidxi = @MEMddbetadxidxi_' filename '_' num2str(S(s)) ';']);
            eval(['Model.exp{s}.dddeltadxidxi = @MEMdddeltadxidxi_' filename '_' num2str(S(s)) ';']);
        end
        eval(['Model.exp{s}.sigma_noise = @MEMsigma_noise_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dsigma_noisedphi = @MEMds_ndp_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.sigma_mean = @MEMsigma_mean_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dsigma_meandphi = @MEMds_mdp_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.sigma_cov = @MEMsigma_cov_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dsigma_covdphi = @MEMds_cdp_' filename '_' num2str(S(s)) ';']);
        if(~onlyPA)
            eval(['Model.exp{s}.ddsigma_noisedphidphi = @MEMdds_ndpdp_' filename '_' num2str(S(s)) ';']);
            eval(['Model.exp{s}.dddsigma_noisedphidphidphi = @MEMddds_ndpdpdp_' filename '_' num2str(S(s)) ';']);
        end
        
        if(Model.integration)
            eval(['Model.exp{s}.ddddsigma_noisedphidphidphidphi = @MEMdddds_ndpdpdpdp_' filename '_' num2str(S(s)) ';']);
        end
        
        eval(['Model.exp{s}.sigma_time = @MEMsigma_time_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dsigma_timedphi = @MEMds_tdp_' filename '_' num2str(S(s)) ';']);
        if(~onlyPA)
            eval(['Model.exp{s}.ddsigma_timedphidphi = @MEMdds_tdpdp_' filename '_' num2str(S(s)) ';']);
            eval(['Model.exp{s}.dddsigma_timedphidphidphi = @MEMddds_tdpdpdp_' filename '_' num2str(S(s)) ';']);
        end
        
        if(Model.integration)
            eval(['Model.exp{s}.ddddsigma_timedphidphidphidphi = @MEMdddds_tdpdpdpdp_' filename '_' num2str(S(s)) ';']);
        end
        
        eval(['Model.exp{s}.phi = @MEMphi_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dphidbeta = @MEMdphidbeta_' filename '_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dphidb = @MEMdphidb_' filename '_' num2str(S(s)) ';']);
        if(~onlyPA)
            eval(['Model.exp{s}.ddphidbdb = @MEMddphidbdb_' filename '_' num2str(S(s)) ';']);
            eval(['Model.exp{s}.ddphidbetadbeta = @MEMddphidbetadbeta_' filename '_' num2str(S(s)) ';']);
            eval(['Model.exp{s}.ddphidbdbeta = @MEMddphidbdbeta_' filename '_' num2str(S(s)) ';']);
        end
    end
end
end
