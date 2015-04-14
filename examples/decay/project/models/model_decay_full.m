function [Model,parameters] = model_decay_full(logbeta,logD,logsigma)
    
    % Parameter struct (for optimization)
    parameters.min = [logbeta;logD;logsigma]-6;
    parameters.max = [logbeta;logD;logsigma]+4;
    parameters.name = {'\mu_{log(m0)}';'\mu_{log(delta)}';...
                   'log(sigma_{m0})';'log(sigma_{m0,delta})';'log(sigma_{delta})'};
    parameters.number = length(parameters.max);
    
    % Covariance parametrisation
    Model.type_D = 'matrix-logarithm';
    
    % Effects Mixture
    Model.param = {'m0','\delta'};
    Model.fixed_effect = [1,1];
    Model.random_effect = [1,1];
    
    % Generate Model
    Model = make_model(Model);
    
end

