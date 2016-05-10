%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_noise(Y,Ym,Sigma,ind)
    Ym = Ym(ind,:);
    Sigma = Sigma(ind,:);
    
if nargout >=1
    % J_D
    varargout{1} = nansum(0.5*((Y - Ym(:))./Sigma(:)).^2 + 0.5*log(2*pi*Sigma(:).^2));
    if nargout >= 3
        % dJ_DdY
        varargout{2} = transpose((Y - Ym(:))./(Sigma(:).^2));
        % dJ_DdSigma
        varargout{3} = transpose(- (((Y - Ym(:)).^2)./(Sigma(:).^3)) + 1./Sigma(:));
        if nargout >= 4
            %ddJ_DdYdY
            varargout{4} = transpose(1./(Sigma(:).^2));
            %ddJ_DdYdSigma
            varargout{5} = transpose(-2*(Y - Ym(:))./(Sigma(:).^3));
            %ddJ_DdSigmadSigma
            varargout{6} = transpose(3*(((Y - Ym(:)).^2)./(Sigma(:).^4)) - 1./(Sigma(:).^2));
            if nargout >= 7
                %dddJ_DdYdYdY
                varargout{7} = transpose(zeros(size(Y)));
                %dddJ_DdYdYdSigma
                varargout{8} = transpose(- 2./(Sigma(:).^3));
                %dddJ_DdYdSigmadSigma
                varargout{9} = transpose(6*(Y - Ym(:))./(Sigma(:).^4));
                %dddJ_DdSigmadSigmadSigma
                varargout{10} = transpose(- 12*(((Y - Ym(:)).^2)./(Sigma(:).^5)) + 2./(Sigma(:).^3));
                if nargout >= 11
                    %ddddJ_DdYdYdYdY
                    varargout{11} = transpose(zeros(size(Y)));
                    %ddddJ_DdYdYdYdSigma
                    varargout{12} = transpose(zeros(size(Y)));
                    %ddddJ_DdYdYdSigmadSigma
                    varargout{13} = transpose(6./(Sigma(:).^4));
                    %ddddJ_DdYdSigmadSigmadSigma
                    varargout{14} = transpose(- 24*((Y - Ym(:))./(Sigma(:).^5)));
                    %ddddJ_DdSigmadSigmadSigmadSigma
                    varargout{15} = transpose(60*(((Y - Ym(:)).^2)./(Sigma(:).^6)) - 6./(Sigma(:).^4));
                end
            end
        end
    end
end
end