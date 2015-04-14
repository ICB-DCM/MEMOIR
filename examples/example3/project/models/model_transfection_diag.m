function [Model,parameters] = model_transfection_diag(logbeta,logD,logsigma)
    
    % Parameter struct (for optimization)
    parameters.min = [logbeta;logD]-6;
    parameters.max = [logbeta;logD]+3;
    parameters.name = {'\beta_{t0}';'\beta_{m_0*k_{TL}}';'\beta_{\beta}';'\beta_{\delta}';...
        '\delta_{11}';'\delta_{22}';'\delta_{33}';'\delta_{44}'};
    parameters.number = length(parameters.max);
    
    % Covariance parametrisation
    Model.type_D = 'diag-matrix-logarithm';
    
    % Effects Mixture
    Model.param = {'t0','m0_kTL','beta','delta'};
    Model.common_effect = [1,1,1,1];
    Model.random_effect = [1,1,1,0];
    
    % Generate Model
    Model = make_model(Model);

end