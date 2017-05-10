function varargout = simulateForSP(model, tout, phi, kappa, scale)

    % Simulate model
    if(nargout<2)
        options_simu.sensi = 0;
        sol = model(tout,phi,kappa,options_simu);
    else
        options_simu.sensi = 1;
        sol = model(tout,phi,kappa,options_simu);
    end

    % Switch for the scale of the simulation
    switch scale
        case 'log'
            varargout{1} = log(sol.y + 1e-16);
            if (nargout>1)
                varargout{2} = bsxfun(@times, sol.sy, 1./(sol.y + 1e-16));
            end
            
        case 'log10'
            varargout{1} = log10(sol.y + 1e-16);
            if (nargout>1)
                varargout{2} = bsxfun(@times, sol.sy, 1 ./ (log(10)*(sol.y + 1e-16)));
            end
            
        case 'lin'
            varargout{1} = sol.y;
            if (nargout>1)
                varargout{2} = sol.sy;
            end
            
    end

    if sol.status < 0
        error('Integration error in simulateForSP.m.')
    end

end