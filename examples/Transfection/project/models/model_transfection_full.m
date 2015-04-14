function [Model,parameters] = model_transfection_red_full(logbeta,logD)
    
    
    % Parameter struct (for optimization)
    parameters.min = [logbeta;logD]-6;
    parameters.max = [logbeta;logD]+3;
    parameters.name = {'\mu_{log(t0)}';'\mu_{log(m_0*k_{TL})}';'\mu_{log(\beta)}';'\mu_{log(\delta)}';...
        'log(D_{11})';'log(D_{12})';'log(D_{22})';'log(D_{31})';'log(D_{32})';...
        'log(D_{33})';'log(D_{41})';'log(D_{42})';'log(D_{43})';'log(D_{44})'};
    parameters.number = length(parameters.max);
    
    % Covariance parametrisation
    Model.type_D = 'matrix-logarithm';
    
    % Effects Mixture
    Model.param = {'t0','m0_kTL','beta','delta'};
    Model.fixed_effect = [1,1,1,1];
    Model.random_effect = [1,1,1,1];
    
    % Generate Model
    Model = make_model(Model);
    
end