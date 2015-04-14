% complete_model uses the symbolic definition of the  parametrisation in 
% the input model struct to generate m-files which can be used to evaluate 
% said symbolic expressions
% 
% USAGE:
% ======
% MODEL = complete_model(MODEL)
%
% INPUTS:
% =======
% Model ... model struct encapsulating the model definition for a MEM
%   .sym ... contains symbolic definition of the overall model
%       .xi ... are the parameter wich are optimised, this usually consist
%       of common effects, the parametrisation of the random effects
%       covariance matrix and the parametrisation of the noise parameters
%       .phi ... is are the mixed effect parametrisation as function of
%       common effects beta and random effects b
%       .beta ...  is the parametrisation of common effects as function of
%       xi
%       .b ... is the parametrisation of random effects
%       .delta ... is the parametrisation of the covariance matrix. this
%       definition should be chosen in accordance to the definition of the
%       respective parametrisation given in Model.type_D
% problem
%
% Outputs:
% ========
% Model ... model struct encapsulating the model definition for a MEM
% problem
%   .exp{s} ... contains all information with respect to experiment number
%   s
%       .N ... number of single cells measured in the experiment
%       .sigma_noise ... noise level in this experiment (used in data
%       generation)
%       .sigma_time ... noise level for event data (used in data
%       generation)
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
filename = [Model.name '_' Model.type_D];

% generate path
mdir  = strrep(mfilename('fullpath'),mfilename('full'),'');

%check for existence of directory
if(~exist([mdir 'MEMfn/' filename ],'dir'))
    mkdir([mdir 'MEMfn/' filename ]);
end

try 
    % check wheter the saved symbolic definition agrees with the current
    % one
    load([mdir 'MEMfn/' filename '/syms.mat'])
    f_xi = isequaln(Model.sym.xi,syms.xi);
    f_phi = isequaln(Model.sym.phi,syms.phi);
    f_beta = isequaln(Model.sym.beta,syms.beta);
    f_b = isequaln(Model.sym.b,syms.b);
    f_delta = isequaln(Model.sym.delta,syms.delta);
    f_exps = length(Model.exp)==length(expsyms);
    for s = 1:length(Model.exp)
        f_phiexp(s) = isequaln(Model.exp{s}.sym.phi,expsyms{s}.phi);
        f_sigma_noiseexp(s) = isequaln(Model.exp{s}.sym.sigma_noise,expsyms{s}.sigma_noise);
        f_sigma_timeexp(s) = isequaln(Model.exp{s}.sym.sigma_time,expsyms{s}.sigma_time);
    end
    if(all([f_xi,f_phi,f_beta,f_b,f_delta,f_exps,f_phiexp,f_sigma_noiseexp,f_sigma_timeexp]))
        if(Model.integration)
            % TBD
        else
            % load old definition
            loadold = true; 
            disp(['Loading previous model definition files!'])
            disp(['To regenerate model, abort and delete ' mdir 'MEMfn/' filename ]);
        end
    end
catch
    disp(['Generating new model definition files!'])
end    

if(~loadold)
    % if we cannot load the old definition, we have to generate a new one

    % remove all other models from the path
    for s=1:length(S)
        if(~strcmp(which(['MEMbeta_' num2str(S(s))]),''))
            if(ispc)
                rmpath(strrep(which(['MEMbeta_' num2str(S(s))]),['\MEMbeta_' num2str(S(s)) '.m'],''));
            else
                rmpath(strrep(which(['MEMbeta_' num2str(S(s))]),['/MEMbeta_' num2str(S(s)) '.m'],''));
            end
        end
    end
    % add the new path
    addpath([mdir 'MEMfn/' filename ]);
    
    % save the symbolic definition as future reference
    syms = Model.sym;
    for s = 1:length(Model.exp)
        expsyms{s} = Model.exp{s}.sym;
    end
    save([mdir 'MEMfn/' filename '/syms.mat'],'syms','expsyms');
    
    % compute number of elements of xi and b
    n_xi = length(Model.sym.xi);
    n_b = length(Model.sym.b);
    
    % construct variance matrix parametrisation
    C = sym('C',[n_b,n_b]);
    switch(Model.type_D)
        case 'matrix-logarithm'
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

        % construct indices for reduced parameters
        Model.exp{s}.ind_beta = find(ismember(Model.sym.beta,symvar(Model.exp{s}.sym.phi)));
        Model.exp{s}.ind_b = find(ismember(Model.sym.b,symvar(Model.exp{s}.sym.phi)));
        Cs = C(Model.exp{s}.ind_b,Model.exp{s}.ind_b);
        Model.exp{s}.ind_delta =  find(ismember(Model.sym.delta,symvar(Cs)));
        
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
        mfun(Model.exp{s}.sym.beta,'file',[mdir 'MEMfn/' filename '/MEMbeta_' num2str(S(s))],'vars',{Model.sym.xi});
        eval(['Model.exp{s}.beta = @MEMbeta_' num2str(S(s)) ';']);
        mfun(Model.exp{s}.sym.delta,'file',[mdir 'MEMfn/' filename '/MEMdelta_' num2str(S(s))],'vars',{Model.sym.xi});
        eval(['Model.exp{s}.delta = @MEMdelta_' num2str(S(s)) ';']);
        
        % dbetadxi
        Model.exp{s}.sym.dbetadxi = simplify(jacobian(Model.exp{s}.sym.beta,Model.sym.xi));
        mfun(Model.exp{s}.sym.dbetadxi,'file',[mdir 'MEMfn/' filename '/MEMdbetadxi_' num2str(S(s))],'vars',{Model.sym.xi});
        eval(['Model.exp{s}.dbetadxi = @MEMdbetadxi_' num2str(S(s)) ';']);
        
        % ddeltadxi
        Model.exp{s}.sym.ddeltadxi = simplify(jacobian(Model.exp{s}.sym.delta,Model.sym.xi));
        mfun(Model.exp{s}.sym.ddeltadxi,'file',[mdir 'MEMfn/' filename '/MEMddeltadxi_' num2str(S(s))],'vars',{Model.sym.xi});
        eval(['Model.exp{s}.ddeltadxi = @MEMddeltadxi_' num2str(S(s)) ';']);
        
        % ddbetadxidxi
        Model.exp{s}.sym.ddbetadxidxi = sym(zeros(n_beta,n_xi,n_xi));
        for j = 1:n_beta
            Model.exp{s}.sym.ddbetadxidxi(j,:,:) = simplify(hessian(Model.exp{s}.sym.beta(j),Model.sym.xi));
        end
        mfun(Model.exp{s}.sym.ddbetadxidxi,'file',[mdir 'MEMfn/' filename '/MEMddbetadxidxi_' num2str(S(s))],'vars',{Model.sym.xi});
        eval(['Model.exp{s}.ddbetadxidxi = @MEMddbetadxidxi_' num2str(S(s)) ';']);
        
        % ddeltadxidxi
        Model.exp{s}.sym.dddeltadxidxi = sym(zeros(n_delta,n_xi,n_xi));
        for j = 1:n_delta
            Model.exp{s}.sym.dddeltadxidxi(j,:,:) = simplify(hessian(Model.exp{s}.sym.delta(j),Model.sym.xi));
        end
        mfun(Model.exp{s}.sym.dddeltadxidxi,'file',[mdir 'MEMfn/' filename '/MEMdddeltadxidxi_' num2str(S(s))],'vars',{Model.sym.xi});
        eval(['Model.exp{s}.dddeltadxidxi = @MEMdddeltadxidxi_' num2str(S(s)) ';']);
        
        
        % sigma_noise(phi)
        mfun(Model.exp{s}.sym.sigma_noise,'file',[mdir 'MEMfn/' filename '/MEMsigma_noise_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.sigma_noise = @MEMsigma_noise_' num2str(S(s)) ';']);
        
        % dsigma_noisedphi
        Model.exp{s}.sym.dsigma_noisedphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
            for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
                Model.exp{s}.sym.dsigma_noisedphi(j,k,:) = jacobian(Model.exp{s}.sym.sigma_noise(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.dsigma_noisedphi,'file',[mdir 'MEMfn/' filename '/MEMdsigma_noisedphi_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.dsigma_noisedphi = @MEMdsigma_noisedphi_' num2str(S(s)) ';']);
        
        % ddsigma_noisedphidphi
        Model.exp{s}.sym.ddsigma_noisedphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi,n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
            for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
                Model.exp{s}.sym.ddsigma_noisedphidphi(j,k,:,:) = hessian(Model.exp{s}.sym.sigma_noise(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.ddsigma_noisedphidphi,'file',[mdir 'MEMfn/' filename '/MEMddsigma_noisedphidphi_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.ddsigma_noisedphidphi = @MEMddsigma_noisedphidphi_' num2str(S(s)) ';']);
        
        % dddsigma_noisedphidphidphi
        Model.exp{s}.sym.dddsigma_noisedphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi,n_phi,n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
            for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
                for m = 1:n_phi
                    Model.exp{s}.sym.dddsigma_noisedphidphidphi(j,k,:,:,m) = diff(Model.exp{s}.sym.ddsigma_noisedphidphi(j,k,:,:),phi(m));
                end
            end
        end
        mfun(Model.exp{s}.sym.dddsigma_noisedphidphidphi,'file',[mdir 'MEMfn/' filename '/MEMdddsigma_noisedphidphidphi_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.dddsigma_noisedphidphidphi = @MEMdddsigma_noisedphidphidphi_' num2str(S(s)) ';']);
        
        if(Model.integration)
            % ddddsigma_noisedphidphidphidphi --- currently disabled due to
            % high computational cost
            Model.exp{s}.sym.ddddsigma_noisedphidphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_noise,1),size(Model.exp{s}.sym.sigma_noise,2),n_phi,n_phi,n_phi,n_phi));
%             for j = 1:size(Model.exp{s}.sym.sigma_noise,1)
%                 for k = 1:size(Model.exp{s}.sym.sigma_noise,2)
%                     for m = 1:n_phi
%                         Model.exp{s}.sym.ddddsigma_noisedphidphidphidphi(j,k,:,:,:,m) = diff(Model.exp{s}.sym.dddsigma_noisedphidphidphi(j,k,:,:,:),phi(m));
%                     end
%                 end
%             end
            mfun(Model.exp{s}.sym.ddddsigma_noisedphidphidphidphi,'file',[mdir 'MEMfn/' filename '/MEMddddsigma_noisedphidphidphidphi_' num2str(S(s))],'vars',{phi});
            eval(['Model.exp{s}.ddddsigma_noisedphidphidphidphi = @MEMddddsigma_noisedphidphidphidphi_' num2str(S(s)) ';']);
        end
        
        % sigma_time(phi)
        if(~isfield(Model.exp{s}.sym,'sigma_time'))
            Model.exp{s}.sym.sigma_time = sym.empty(0,1);
        end
        mfun(Model.exp{s}.sym.sigma_time,'file',[mdir 'MEMfn/' filename '/MEMsigma_time_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.sigma_time = @MEMsigma_time_' num2str(S(s)) ';']);
        
        % dsigma_timedphi
        Model.exp{s}.sym.dsigma_timedphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_time,1)
            for k = 1:size(Model.exp{s}.sym.sigma_time,2)
                Model.exp{s}.sym.dsigma_timedphi(j,k,:) = jacobian(Model.exp{s}.sym.sigma_time(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.dsigma_timedphi,'file',[mdir 'MEMfn/' filename '/MEMdsigma_timedphi_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.dsigma_timedphi = @MEMdsigma_timedphi_' num2str(S(s)) ';']);
        
        % ddsigma_timedphidphi
        Model.exp{s}.sym.ddsigma_timedphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi,n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_time,1)
            for k = 1:size(Model.exp{s}.sym.sigma_time,2)
                Model.exp{s}.sym.ddsigma_timedphidphi(j,k,:,:) = hessian(Model.exp{s}.sym.sigma_time(j,k),phi);
            end
        end
        mfun(Model.exp{s}.sym.ddsigma_timedphidphi,'file',[mdir 'MEMfn/' filename '/MEMddsigma_timedphidphi_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.ddsigma_timedphidphi = @MEMddsigma_timedphidphi_' num2str(S(s)) ';']);
        
        % dddsigma_timedphidphidphi
        Model.exp{s}.sym.dddsigma_timedphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi,n_phi,n_phi));
        for j = 1:size(Model.exp{s}.sym.sigma_time,1)
            for k = 1:size(Model.exp{s}.sym.sigma_time,2)
                for m = 1:n_phi
                    Model.exp{s}.sym.dddsigma_timedphidphidphi(j,k,:,:,m) = diff(Model.exp{s}.sym.ddsigma_timedphidphi(j,k,:,:),phi(m));
                end
            end
        end
        mfun(Model.exp{s}.sym.dddsigma_timedphidphidphi,'file',[mdir 'MEMfn/' filename '/MEMdddsigma_timedphidphidphi_' num2str(S(s))],'vars',{phi});
        eval(['Model.exp{s}.dddsigma_timedphidphidphi = @MEMdddsigma_timedphidphidphi_' num2str(S(s)) ';']);
        
        if(Model.integration)
            % ddddsigma_timedphidphidphidphi --- currently disabled due to
            % high computational cost
            Model.exp{s}.sym.ddddsigma_timedphidphidphidphi = sym(zeros(size(Model.exp{s}.sym.sigma_time,1),size(Model.exp{s}.sym.sigma_time,2),n_phi,n_phi,n_phi,n_phi));
%             for j = 1:size(Model.exp{s}.sym.sigma_time,1)
%                 for k = 1:size(Model.exp{s}.sym.sigma_time,2)
%                     for m = 1:n_phi
%                         Model.exp{s}.sym.ddddsigma_timedphidphidphidphi(j,k,:,:,:,m) = diff(Model.exp{s}.sym.dddsigma_timedphidphidphi(j,k,:,:,:),phi(m));
%                     end
%                 end
%             end
            mfun(Model.exp{s}.sym.ddddsigma_timedphidphidphidphi,'file',[mdir 'MEMfn/' filename '/MEMddddsigma_timedphidphidphidphi_' num2str(S(s))],'vars',{phi});
            eval(['Model.exp{s}.ddddsigma_timedphidphidphidphi = @MEMddddsigma_timedphidphidphidphi_' num2str(S(s)) ';']);
        end
        
        % phi
        mfun(Model.exp{s}.sym.phi,'file',[mdir 'MEMfn/' filename '/MEMphi_' num2str(S(s))],'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.phi = @MEMphi_' num2str(S(s)) ';']);
        
        % dphidbeta
        Model.exp{s}.sym.dphidbeta = simplify(jacobian(Model.exp{s}.sym.phi,Model.exp{s}.sym.beta));
        mfun(Model.exp{s}.sym.dphidbeta,'file',[mdir 'MEMfn/' filename '/MEMdphidbeta_' num2str(S(s))],'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.dphidbeta = @MEMdphidbeta_' num2str(S(s)) ';']);
        
        % dphidb
        Model.exp{s}.sym.dphidb = simplify(jacobian(Model.exp{s}.sym.phi,Model.exp{s}.sym.b));
        mfun(Model.exp{s}.sym.dphidb,'file',[mdir 'MEMfn/' filename '/MEMdphidb_' num2str(S(s))],'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.dphidb = @MEMdphidb_' num2str(S(s)) ';']);
        
        
        % ddphidbdb
        Model.exp{s}.sym.ddphidbdb  = sym(zeros(n_phi,n_b,n_b));
        for j = 1:n_phi
            Model.exp{s}.sym.ddphidbdb(j,:,:) = simplify(hessian(Model.exp{s}.sym.phi(j),Model.exp{s}.sym.b));
        end
        mfun(Model.exp{s}.sym.ddphidbdb,'file',[mdir 'MEMfn/' filename '/MEMddphidbdb_' num2str(S(s))],'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.ddphidbdb = @MEMddphidbdb_' num2str(S(s)) ';']);
        
        % dphidbetadbeta
        Model.exp{s}.sym.ddphidbetadbeta  = sym(zeros(n_phi,n_beta,n_beta));
        for j = 1:n_phi
            Model.exp{s}.sym.ddphidbetadbeta(j,:,:) = simplify(hessian(Model.exp{s}.sym.phi(j),Model.exp{s}.sym.beta));
        end
        mfun(Model.exp{s}.sym.ddphidbetadbeta,'file',[mdir 'MEMfn/' filename '/MEMddphidbetadbeta_' num2str(S(s))],'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.ddphidbetadbeta = @MEMddphidbetadbeta_' num2str(S(s)) ';']);
        
        % dphidbdbeta
        Model.exp{s}.sym.ddphidbetadb  = sym(zeros(n_phi,n_beta,n_b));
        for j = 1:n_phi
            Model.exp{s}.sym.ddphidbdbeta(j,:,:) = simplify(jacobian(jacobian(Model.exp{s}.sym.phi(j),Model.exp{s}.sym.b),Model.exp{s}.sym.beta));
        end
        mfun(Model.exp{s}.sym.ddphidbdbeta,'file',[mdir 'MEMfn/' filename '/MEMddphidbdbeta_' num2str(S(s))],'vars',{Model.exp{s}.sym.beta,Model.exp{s}.sym.b});
        eval(['Model.exp{s}.ddphidbdbeta = @MEMddphidbdbeta_' num2str(S(s)) ';']);
        
    end

else
    % if we can load the old definition, we just have to attach the m-files
    % to the model struct
    
    % remove all other models from the path
    for s=1:length(S)
        if(~strcmp(which(['MEMbeta_' num2str(S(s))]),''))
            if(ispc)
                rmpath(strrep(which(['MEMbeta_' num2str(S(s))]),['\MEMbeta_' num2str(S(s)) '.m'],''));
            else
                rmpath(strrep(which(['MEMbeta_' num2str(S(s))]),['/MEMbeta_' num2str(S(s)) '.m'],''));
            end
        end
    end
    % add the new path
    addpath([mdir 'MEMfn/' filename ]);
    
    % save the symbolic definition as future reference
    syms = load([mdir 'MEMfn/' filename '/syms.mat']);
    Model.sym = syms.syms;
    
    % compute number of elements of xi and b
    n_xi = length(Model.sym.xi);
    n_b = length(Model.sym.b);
    
    
    % construct variance matrix parametrisation
    C = sym('C',[n_b,n_b]);
    switch(Model.type_D)
        case 'matrix-logarithm'
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
        Model.exp{s}.ind_beta = find(ismember(Model.sym.beta,symvar(Model.exp{s}.sym.phi)));
        Model.exp{s}.ind_b = find(ismember(Model.sym.b,symvar(Model.exp{s}.sym.phi)));
        Cs = C(Model.exp{s}.ind_b,Model.exp{s}.ind_b);
        Model.exp{s}.ind_delta =  find(ismember(Model.sym.delta,symvar(Cs)));
        
        % constructe reduced parameters
        Model.exp{s}.sym.beta = Model.sym.beta(Model.exp{s}.ind_beta);
        Model.exp{s}.sym.b = Model.sym.b(Model.exp{s}.ind_b);
        Model.exp{s}.sym.delta = Model.sym.delta(Model.exp{s}.ind_delta);
        
        eval(['Model.exp{s}.beta = @MEMbeta_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.delta = @MEMdelta_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dbetadxi = @MEMdbetadxi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddeltadxi = @MEMddeltadxi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddbetadxidxi = @MEMddbetadxidxi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dddeltadxidxi = @MEMdddeltadxidxi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.sigma_noise = @MEMsigma_noise_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dsigma_noisedphi = @MEMdsigma_noisedphi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddsigma_noisedphidphi = @MEMddsigma_noisedphidphi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dddsigma_noisedphidphidphi = @MEMdddsigma_noisedphidphidphi_' num2str(S(s)) ';']);
        
        if(Model.integration)
            eval(['Model.exp{s}.ddddsigma_noisedphidphidphidphi = @MEMddddsigma_noisedphidphidphidphi_' num2str(S(s)) ';']);
        end
        
        eval(['Model.exp{s}.sigma_time = @MEMsigma_time_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dsigma_timedphi = @MEMdsigma_timedphi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddsigma_timedphidphi = @MEMddsigma_timedphidphi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dddsigma_timedphidphidphi = @MEMdddsigma_timedphidphidphi_' num2str(S(s)) ';']);
        
        if(Model.integration)
            eval(['Model.exp{s}.ddddsigma_timedphidphidphidphi = @MEMddddsigma_timedphidphidphidphi_' num2str(S(s)) ';']);
        end
        
        eval(['Model.exp{s}.phi = @MEMphi_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dphidbeta = @MEMdphidbeta_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.dphidb = @MEMdphidb_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddphidbdb = @MEMddphidbdb_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddphidbetadbeta = @MEMddphidbetadbeta_' num2str(S(s)) ';']);
        eval(['Model.exp{s}.ddphidbdbeta = @MEMddphidbdbeta_' num2str(S(s)) ';']);
    end
end
end
