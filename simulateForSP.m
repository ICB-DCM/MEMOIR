function varargout = simulateForSP(model,tout,phi,kappa)

% Simulate model
if(nargout<2)
    options_simu.sensi = 0;
    sol = model(tout,phi,kappa,[],options_simu);
    varargout{1} = log(sol.y);
else
    options_simu.sensi = 1;
    sol = model(tout,phi,kappa,[],options_simu);
    varargout{1} = log(sol.y);
    varargout{2} = bsxfun(@times,sol.sy,1./sol.y);
end

if sol.status < 0
    error('Integration error in simulateForSP.m.')
end


end