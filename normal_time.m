%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_time(T,Tm,R,Sigma,ind)
    
    
if nargout >=1
    % J_T
    varargout{1} = sum(0.5*((T(ind) - Tm(ind))./Sigma(ind)).^2 + 0.5*((R(ind))./Sigma(ind)).^2 + log(2*pi*Sigma(ind).^2));
    if nargout >= 2
        % dJ_TdT
        varargout{2} = transpose((T(ind) - Tm(ind))./(Sigma(ind).^2));
        % dJ_TdR
        varargout{3} = transpose(R(ind)./(Sigma(ind).^2));
        % dJ_TdSigma
        varargout{4} = transpose(- (((T(ind) - Tm(ind)).^2)./(Sigma(ind).^3)) - (((R(ind)).^2)./(Sigma(ind).^3)) + 2./Sigma(ind));
        if nargout >= 5
            %ddJ_TdTdT
            varargout{5} = transpose(1./(Sigma(ind).^2));
            %ddJ_TdTdR
            varargout{6} = transpose(zeros(size(T(ind))));
            %ddJ_TdRdR
            varargout{7} = transpose(1./(Sigma(ind).^2));
            %ddJ_TdTdSigma
            varargout{8} = transpose(-2*(T(ind) - Tm(ind))./(Sigma(ind).^3));
            %ddJ_TdRdSigma
            varargout{9} = transpose(-2*(R(ind))./(Sigma(ind).^3));
            %ddJ_TdSigmadSigma
            varargout{10} = transpose(3*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^4)) + 3*(((R(ind)).^2)./(Sigma(ind).^4)) - 2./(Sigma(ind).^2));
            if nargout >= 11
                %dddJ_TdTdTdT
                varargout{11} = transpose(zeros(size(T(ind))));
                %dddJ_TdTdTdR
                varargout{12} = transpose(zeros(size(T(ind))));
                %dddJ_TdTdRdR
                varargout{13} = transpose(zeros(size(T(ind))));
                %dddJ_TdRdRdR
                varargout{14} = transpose(zeros(size(T(ind))));
                %dddJ_TdTdTdSigma
                varargout{15} = transpose(- 2./(Sigma(ind).^3));
                %dddJ_TdTdRdSigma
                varargout{16} = transpose(zeros(size(T(ind))));
                %dddJ_TdRdRdSigma
                varargout{17} = transpose(- 2./(Sigma(ind).^3));
                %dddJ_TdTdSigmadSigma
                varargout{18} = transpose(6*(T(ind) - Tm(ind))./(Sigma(ind).^4));
                %dddJ_TdRdSigmadSigma
                varargout{19} = transpose(6*(R(ind))./(Sigma(ind).^4));
                %dddJ_TdSigmadSigmadSigma
                varargout{20} = transpose(- 12*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^5)) - 12*(((R(ind)).^2)./(Sigma(ind).^5)) + 4./(Sigma(ind).^3));
                if nargout >= 20
                    %ddddJ_TdTdTdTdT
                    varargout{21} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdTdTdR
                    varargout{22} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdTdRdR
                    varargout{23} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdRdRdR
                    varargout{24} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdRdRdRdR
                    varargout{25} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdTdTdSigma
                    varargout{26} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdTdRdSigma
                    varargout{27} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdRdRdSigma
                    varargout{28} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdRdRdRdSigma
                    varargout{29} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdTdTdSigmadSigma
                    varargout{30} = transpose(6./(Sigma(ind).^4));
                    %ddddJ_TdTdRdSigmadSigma
                    varargout{31} = transpose(zeros(size(T(ind))));
                    %ddddJ_TdRdRdSigmadSigma
                    varargout{32} = transpose(6./(Sigma(ind).^4));
                    %ddddJ_TdTdSigmadSigmadSigma
                    varargout{33} = transpose(- 24*((T(ind) - Tm(ind))./(Sigma(ind).^5)));
                    %ddddJ_TdRdSigmadSigmadSigma
                    varargout{34} = transpose(- 24*((R(ind))./(Sigma(ind).^5)));
                    %ddddJ_TdSigmadSigmadSigmadSigma
                    varargout{35} = transpose(60*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^6)) + 60*(((R(ind)).^2)./(Sigma(ind).^6)) - 12./(Sigma(ind).^4));
                end
            end
        end
    end
end
end