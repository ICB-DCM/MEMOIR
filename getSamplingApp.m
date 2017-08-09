% Function of the toolbox MEMOIR
% getSampleApp.m computes the integral in mixed effect models by sampling
% according to some initial distribution, of which mean and covariance
% matrix were passed to it.

function SamplingStruct = getSamplingApp(nonlinFun, xi, ExpStruct, optionsSampling)
    
    % Check how many derivatives are to be computed
    nDeriv = optionsSampling.nDeriv;
    if (nDeriv > 0)
        error('Currently, derivatives are not supported within the sampling routine!');
    end
    
    % Preallocate the simulation results
    nSamples = optionsSampling.nSamples;
    Y = nan(length(ExpStruct.t), 1, nSamples);
    
    % Get parameter vectors for the current experiment from xi
    beta = ExpStruct.beta(xi);
    delta = ExpStruct.delta(xi);

    % Compute covariance matrix for t = 0
    [D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,'diag-matrix-logarithm');
    n_b = size(D,1);
    
    % Create sampels
    samples = transpose(mvnrnd(zeros(n_b,1), D, nSamples));
    
    % propagate sampples through model
    for iSample = 1 : nSamples
        Y(:, :, iSample) = nonlinFun(ExpStruct.phi(beta,samples(:,iSample)));
    end
    
    my = mean(Y,3);
    
    SamplingStruct = struct(...
        'Y', Y, ...
        'my', my);

end