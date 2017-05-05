function varargout = simulateForSP(model,tout,phi,kappa)

% Simulate model
if(nargout<2)
    options_simu.sensi = 0;
    sol = model(tout,phi,kappa,options_simu);
    varargout{1} = sol.y; % sol.y+1e-16;
else
    options_simu.sensi = 1;
    sol = model(tout,phi,kappa,options_simu);
    varargout{1} = sol.y; % log(sol.y + 1e-16);
    varargout{2} = sol.sy; % bsxfun(@times, sol.sy, 1./(sol.y + 1e-16));
end

% if(nargout<2)
%     options_simu.sensi = 0;
%     sol = model(tout,phi,kappa,options_simu);
%     varargout{1} = sol.y;
% else
%     options_simu.sensi = 1;
%     sol = model(tout,phi,kappa,options_simu);
%     varargout{1} = sol.y;
%     varargout{2} = sol.sy;
% end

if sol.status < 0
    error('Integration error in simulateForSP.m.')
end


end