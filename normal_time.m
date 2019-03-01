%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_T = normal_time(T,Tm,R,Sigma,ind,nderiv)

T = T(ind);
Tm = Tm(ind);
R = R(ind);
Sigma = Sigma(ind);

if nderiv >= 0
    % J_T
    J_T.val = sum(0.5*((T - Tm)./Sigma).^2 + 0.5*((R)./Sigma).^2 + log(2*pi*Sigma.^2));
    if nderiv >= 1
        % dJ_TdT
        J_T.dT = transpose((T - Tm)./(Sigma.^2));
        % dJ_TdR 
            J_T.dR = transpose(R./(Sigma.^2));
% % %         end
        % dJ_TdSigma
        J_T.dSigma = transpose(- (((T - Tm).^2)./(Sigma.^3)) - (((R).^2)./(Sigma.^3)) + 2./Sigma);
        if nderiv >= 2
            %ddJ_TdTdT
            J_T.dTdT = diag(1./(Sigma.^2));
            %ddJ_TdTdR
            J_T.dTdR = diag(zeros(size(T)));
            %ddJ_TdRdR
            J_T.dRdR = diag(1./(Sigma.^2));
            %ddJ_TdTdSigma
            J_T.dTdSigma = diag(-2*(T - Tm)./(Sigma.^3));
            %ddJ_TdRdSigma
            J_T.dRdSigma = diag(-2*(R)./(Sigma.^3));
            %ddJ_TdSigmadSigma
            J_T.dSigmadSigma = diag(3*(((T - Tm).^2)./(Sigma.^4)) + 3*(((R).^2)./(Sigma.^4)) - 2./(Sigma.^2));
            if nderiv >= 3
                %dddJ_TdTdTdT
                J_T.dTdTdT = transpose(zeros(size(T)));
                %dddJ_TdTdTdR
                J_T.dTdTdR = transpose(zeros(size(T)));
                %dddJ_TdTdRdR
                J_T.dTdRdR = transpose(zeros(size(T)));
                %dddJ_TdRdRdR
                J_T.dRdRdR = transpose(zeros(size(T)));
                %dddJ_TdTdTdSigma
                J_T.dTdTdSigma = transpose(- 2./(Sigma.^3));
                %dddJ_TdTdRdSigma
                J_T.dTdRdSigma = transpose(zeros(size(T)));
                %dddJ_TdRdRdSigma
                J_T.dRdRdSigma = transpose(- 2./(Sigma.^3));
                %dddJ_TdTdSigmadSigma
                J_T.dTdSigmadSigma = transpose(6*(T - Tm)./(Sigma.^4));
                %dddJ_TdRdSigmadSigma
                J_T.dRdSigmadSigma = transpose(6*(R)./(Sigma.^4));
                %dddJ_TdSigmadSigmadSigma
                J_T.dSigmadSigmadSigma = transpose(- 12*(((T - Tm).^2)./(Sigma.^5)) - 12*(((R).^2)./(Sigma.^5)) + 4./(Sigma.^3));
%                 if nargout >= 20
%                     %ddddJ_TdTdTdTdT
%                     varargout{21} = transpose(zeros(size(T)));
%                     if(isempty(varargout{21}))
%                         varargout{21} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdTdTdR
%                     varargout{22} = transpose(zeros(size(T)));
%                     if(isempty(varargout{22}))
%                         varargout{22} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdTdRdR
%                     varargout{23} = transpose(zeros(size(T)));
%                     if(isempty(varargout{23}))
%                         varargout{23} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdRdRdR
%                     varargout{24} = transpose(zeros(size(T)));
%                     if(isempty(varargout{24}))
%                         varargout{24} = zeros(1,0);
%                     end
%                     %ddddJ_TdRdRdRdR
%                     varargout{25} = transpose(zeros(size(T)));
%                     if(isempty(varargout{25}))
%                         varargout{25} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdTdTdSigma
%                     varargout{26} = transpose(zeros(size(T)));
%                     if(isempty(varargout{26}))
%                         varargout{26} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdTdRdSigma
%                     varargout{27} = transpose(zeros(size(T)));
%                     if(isempty(varargout{27}))
%                         varargout{27} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdRdRdSigma
%                     varargout{28} = transpose(zeros(size(T)));
%                     if(isempty(varargout{28}))
%                         varargout{28} = zeros(1,0);
%                     end
%                     %ddddJ_TdRdRdRdSigma
%                     varargout{29} = transpose(zeros(size(T)));
%                     if(isempty(varargout{29}))
%                         varargout{29} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdTdSigmadSigma
%                     varargout{30} = transpose(6./(Sigma.^4));
%                     if(isempty(varargout{30}))
%                         varargout{30} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdRdSigmadSigma
%                     varargout{31} = transpose(zeros(size(T)));
%                     if(isempty(varargout{31}))
%                         varargout{31} = zeros(1,0);
%                     end
%                     %ddddJ_TdRdRdSigmadSigma
%                     varargout{32} = transpose(6./(Sigma.^4));
%                     if(isempty(varargout{32}))
%                         varargout{32} = zeros(1,0);
%                     end
%                     %ddddJ_TdTdSigmadSigmadSigma
%                     varargout{33} = transpose(- 24*((T - Tm)./(Sigma.^5)));
%                     if(isempty(varargout{33}))
%                         varargout{33} = zeros(1,0);
%                     end
%                     %ddddJ_TdRdSigmadSigmadSigma
%                     varargout{34} = transpose(- 24*((R)./(Sigma.^5)));
%                     if(isempty(varargout{34}))
%                         varargout{34} = zeros(1,0);
%                     end
%                     %ddddJ_TdSigmadSigmadSigmadSigma
%                     varargout{35} = transpose(60*(((T - Tm).^2)./(Sigma.^6)) + 60*(((R).^2)./(Sigma.^6)) - 12./(Sigma.^4));
%                     if(isempty(varargout{35}))
%                         varargout{35} = zeros(1,0);
%                     end
%                 end
            end
        end
    end
end
end