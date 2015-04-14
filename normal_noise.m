%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_noise(Y,Ym,Sigma,ind)
if nargout >=1
    % J_D
    varargout{1} = sum(0.5*((Y(ind) - Ym(ind))./Sigma(ind)).^2 + 0.5*log(sqrt(2*pi)*Sigma(ind).^2));
    if nargout >= 3
        % dJ_DdY
        varargout{2} = transpose((Y(ind) - Ym(ind))./(Sigma(ind).^2));
        % dJ_DdSigma
        varargout{3} = transpose(- (((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^3)) + 1./Sigma(ind));
        if nargout >= 4
            %ddJ_DdYdY
            varargout{4} = transpose(1./(Sigma(ind).^2));
            %ddJ_DdYdSigma
            varargout{5} = transpose(-2*(Y(ind) - Ym(ind))./(Sigma(ind).^3));
            %ddJ_DdSigmadSigma
            varargout{6} = transpose(3*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^4)) - 1./(Sigma(ind).^2));
            if nargout >= 7
                %dddJ_DdYdYdY
                varargout{7} = transpose(zeros(size(Y(ind))));
                %dddJ_DdYdYdSigma
                varargout{8} = transpose(- 2./(Sigma(ind).^3));
                %dddJ_DdYdSigmadSigma
                varargout{9} = transpose(6*(Y(ind) - Ym(ind))./(Sigma(ind).^4));
                %dddJ_DdSigmadSigmadSigma
                varargout{10} = transpose(- 12*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^5)) + 2./(Sigma(ind).^3));
                if nargout >= 11
                    %ddddJ_DdYdYdYdY
                    varargout{11} = transpose(zeros(size(Y(ind))));
                    %ddddJ_DdYdYdYdSigma
                    varargout{12} = transpose(zeros(size(Y(ind))));
                    %ddddJ_DdYdYdSigmadSigma
                    varargout{13} = transpose(6./(Sigma(ind).^4));
                    %ddddJ_DdYdSigmadSigmadSigma
                    varargout{14} = transpose(- 24*((Y(ind) - Ym(ind))./(Sigma(ind).^5)));
                    %ddddJ_DdSigmadSigmadSigmadSigma
                    varargout{15} = transpose(60*(((Y(ind) - Ym(ind)).^2)./(Sigma(ind).^6)) - 6./(Sigma(ind).^4));
                end
            end
        end
    end
end
end