% make_model uses the information contained in the model struct to generate
% symbolic expressions for parameters
%
% USAGE:
% ======
% MODEL = make_model(MODEL)
%
% INPUTS:
% =======
% Model ... model struct encapsulating the model definition for a MEM
% problem
% the model struct must have the following fields of equal length
%   .param cell array of parameter names
%   .common_effect indicates whether there is a common effect
%   component for that parameter. common effects components are shared
%   across individuals
%   .random_effect indicates whether there is a random effect
%   component for that parameter. random effect components are different
%   for every individual
%
% Outputs:
% ========
% Model ... model struct encapsulating the model definition for a MEM
% problem
% the function make_model will add the following fields to the struct
%   .sym ... contains symbolic definition of the overall model
%       .xi ... are the parameter wich are optimised, this usually consist
%       of common effects, the parametrisation of the random effects
%       covariance matrix and the parametrisation of the noise parameters
%       .phi ... is are the mixed effect parametrisation as function of
%       common effects beta and random effects b
%       .beta ...  is the parametrisation of common effects as function of
%       xi
%       .b ... is the parametrisation of random effects
%       .delta ... is the parametrisation of the covariance matrix of the
%       random effect.
% 2015/04/14 Fabian Froehlich

function [ Model ] = make_model( Model )

% backwards compatibility
if(isfield(Model,'fixed_effect') && ~isfield(Model,'common_effect'))
    Model.common_effect = Model.fixed_effect;
end

if(any([length(Model.param)~=length(Model.common_effect),length(Model.param)~=length(Model.random_effect),length(Model.common_effect)~=length(Model.random_effect)]))
    error('Size of Model.param, Model.common_effect and Model.random_effect (Model.fixed_effect) do not agree')
end

% strip parameters of bad characters as this will cause problems with generation of symbolic variables
for j = 1:length(Model.param)
    Model.param{j} = strrep(Model.param{j},' ','_');
    Model.param{j} = strrep(Model.param{j},'\','');
end

switch(Model.type_D)
    case 'matrix-logarithm'
        n = length(Model.param);
        n_f = sum(Model.common_effect);
        n_r = sum(Model.random_effect);
        % initialise symbolic variables
        Model.sym.xi = sym('xi',[n_f+(n_r^2+n_r)/2,1]);
        Model.sym.phi = sym('phi',[sum(arrayfun(@or,Model.common_effect,Model.random_effect)),1]);
        Model.sym.beta = sym('beta',[n_f,1]);
        Model.sym.b = sym('b',[n_r,1]);
        Model.sym.delta = sym('delta',[(n_r^2+n_r)/2,1]);
        k = 0;
        k_f = 0;
        k_r = 0;
        for j = 1:n;
            if(Model.common_effect(j))
                k_f = k_f + 1;
                % construct variables corresponding to common effects
                eval(['syms M_' Model.param{j} ';']);
                eval(['Model.sym.xi(k_f) = M_' Model.param{j} ';']);
                eval(['Model.sym.beta(k_f) = M_' Model.param{j} ';']);
            end
            if(Model.random_effect(j))
                % construct variables corresponding to radom effects
                k_r = k_r + 1;
                eval(['syms b_' Model.param{j} ';']);
                ll = find(Model.random_effect(1:k_r));
                for l = 1:length(ll)
                    m = ll(l);
                    eval(['syms C_' Model.param{j} '_' Model.param{m} ';']);
                    eval(['Model.sym.xi(sum(Model.common_effect)+((k_r-1)^2+(k_r-1))/2+l) = C_' Model.param{j} '_' Model.param{m} ';']);
                    eval(['Model.sym.delta(((k_r-1)^2+(k_r-1))/2+l) = C_' Model.param{j} '_' Model.param{m} ';'])
                end
                eval(['Model.sym.b(k_r) = b_' Model.param{j} ';']);
            end
            if(or(Model.common_effect(j),Model.random_effect(j)))
                % construct variables corresponding to mixed effects
                k = k + 1;
                eval(['syms p_' Model.param{j} ';']);
                eval(['Model.sym.phi(k) = p_' Model.param{j} ';']);
            end
        end
    case 'diag-matrix-logarithm'
        n = length(Model.param);
        n_f = sum(Model.common_effect);
        n_r = sum(Model.random_effect);
        
        % initialise symbolic variables
        Model.sym.xi = sym('xi',[n_f+n_r,1]);
        Model.sym.phi = sym('phi',[sum(arrayfun(@or,Model.common_effect,Model.random_effect)),1]);
        Model.sym.beta = sym('beta',[n_f,1]);
        Model.sym.b = sym('b',[n_r,1]);
        Model.sym.delta = sym('delta',[n_r,1]);
        k = 0;
        k_f = 0;
        k_r = 0;
        for j = 1:n;
            if(Model.common_effect(j))
                % construct variables corresponding to common effects
                k_f = k_f + 1;
                eval(['syms M_' Model.param{j} ';']);
                eval(['Model.sym.xi(k_f) = M_' Model.param{j} ';']);
                eval(['Model.sym.beta(k_f) = M_' Model.param{j} ';']);
                
            end
            if(Model.random_effect(j))
                % construct variables corresponding to radom effects
                k_r = k_r + 1;
                eval(['syms C_' Model.param{j} ' b_' Model.param{j} ';']);
                eval(['Model.sym.xi(sum(Model.common_effect)+k_r) = C_' Model.param{j} ';']);
                eval(['Model.sym.b(k_r) = b_' Model.param{j} ';']);
                eval(['Model.sym.delta(k_r) = C_' Model.param{j} ';']);
            end
            if(or(Model.common_effect(j),Model.random_effect(j)))
                % construct variables corresponding to mixed effects
                k = k + 1;
                eval(['syms p_' Model.param{j} ';']);
                eval(['Model.sym.phi(k) = p_' Model.param{j} ';']);
            end
        end
end
end

