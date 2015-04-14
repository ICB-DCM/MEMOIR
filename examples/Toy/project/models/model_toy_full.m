function [Model,parameters] = model_toy_full(logbeta,logD,logsigma)

    % Parameter struct (for optimization)
    parameters.min = [logbeta;logD;logsigma]-6;
    parameters.max = [logbeta;logD;logsigma]+4;
    parameters.name = {'\mu_{1}';'\mu_{2}';...
                   'log(D_{11})';'log(D_{12}';'log(D_{22})'};
    parameters.number = length(parameters.max);
    
    % Covariance parametrisation
    Model.type_D = 'matrix-logarithm';
    
    % Effects Mixture
    Model.param = {'mu1','mu2'};
    Model.fixed_effect = [1,1];
    Model.random_effect = [1,1];
    
    % Generate Model
    Model = make_model(Model);
end

