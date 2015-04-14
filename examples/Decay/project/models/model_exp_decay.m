function [Model,parameters] = model_exp_decay(logbeta,logD,logsigma)
    
    % Parameter struct (for optimization)
    parameters.min = [logbeta;logD;logsigma]-6;
    parameters.max = [logbeta;logD;logsigma]+4;
    parameters.name = {'\mu_{log(x0)}';'\mu_{log(d)}';...
                   'log(sigma_d)'};
    parameters.number = length(parameters.max);
    
    % Covariance parametrisation
    Model.type_D = 'diag-matrix-logarithm';
    
    % Effects Mixture
    Model.param = {'x0','d'};
    Model.fixed_effect = [1,1];
    Model.random_effect = [0,1];
    
    % Generate Model
    Model = make_model(Model);
    
end

