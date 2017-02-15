function varargout = objective_phi_wrapper(phi,Model,Data,s,i)

nderiv = nargout-1;
FIMflag = nargout >=2;

[J_D,J_T] = objective_phi(Model.exp{s},Data{s},phi,s,i,[],nderiv,FIMflag);

if nargout >= 1;
    varargout{1} = -J_D.val - J_T.val;
    if nargout >= 2;
        varargout{2} = -J_D.dphi - J_T.dphi;
        if nargout >= 3;
            varargout{3} = -J_D.FIM - J_T.FIM;
        end
    end
end




end