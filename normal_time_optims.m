%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR NOISE MODELS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = normal_time(T,Tm,R,Sigma,ind)
    
    varargout = cell(nargout,1); % whyever that is necessary ...
    if nargout >=1
        % J_T
        varargout{1} = sum(0.5*((T(ind) - Tm(ind))./Sigma(ind)).^2 + 0.5*((R(ind))./Sigma(ind)).^2 + log(2*pi*Sigma(ind).^2));
        if nargout >= 2
            % dJ_TdT
            varargout{2} = transpose((T(ind) - Tm(ind))./(Sigma(ind).^2));
            if(isempty(varargout{2}))
                varargout{2} = zeros(1,0);
            end
            % dJ_TdR
            varargout{3} = transpose(R(ind)./(Sigma(ind).^2));
            if(isempty(varargout{3}))
                varargout{3} = zeros(1,0);
            end
            % dJ_TdSigma
            varargout{4} = transpose(- (((T(ind) - Tm(ind)).^2)./(Sigma(ind).^3)) - (((R(ind)).^2)./(Sigma(ind).^3)) + 2./Sigma(ind));
            if(isempty(varargout{4}))
                varargout{4} = zeros(1,0);
            end
            if nargout >= 5
                %ddJ_TdTdT
                varargout{5} = transpose(1./(Sigma(ind).^2));
                if(isempty(varargout{5}))
                    varargout{5} = zeros(1,0);
                end
                %ddJ_TdTdR
                varargout{6} = transpose(zeros(size(T(ind))));
                if(isempty(varargout{6}))
                    varargout{6} = zeros(1,0);
                end
                %ddJ_TdRdR
                varargout{7} = transpose(1./(Sigma(ind).^2));
                if(isempty(varargout{7}))
                    varargout{7} = zeros(1,0);
                end
                %ddJ_TdTdSigma
                varargout{8} = transpose(-2*(T(ind) - Tm(ind))./(Sigma(ind).^3));
                if(isempty(varargout{8}))
                    varargout{8} = zeros(1,0);
                end
                %ddJ_TdRdSigma
                varargout{9} = transpose(-2*(R(ind))./(Sigma(ind).^3));
                if(isempty(varargout{9}))
                    varargout{9} = zeros(1,0);
                end
                %ddJ_TdSigmadSigma
                varargout{10} = transpose(3*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^4)) + 3*(((R(ind)).^2)./(Sigma(ind).^4)) - 2./(Sigma(ind).^2));
                if(isempty(varargout{10}))
                    varargout{10} = zeros(1,0);
                end
                if nargout >= 11
                    %dddJ_TdTdTdT
                    varargout{11} = transpose(zeros(size(T(ind))));
                    if(isempty(varargout{11}))
                        varargout{11} = zeros(1,0);
                    end
                    %dddJ_TdTdTdR
                    varargout{12} = transpose(zeros(size(T(ind))));
                    if(isempty(varargout{12}))
                        varargout{12} = zeros(1,0);
                    end
                    %dddJ_TdTdRdR
                    varargout{13} = transpose(zeros(size(T(ind))));
                    if(isempty(varargout{13}))
                        varargout{13} = zeros(1,0);
                    end
                    %dddJ_TdRdRdR
                    varargout{14} = transpose(zeros(size(T(ind))));
                    if(isempty(varargout{14}))
                        varargout{14} = zeros(1,0);
                    end
                    %dddJ_TdTdTdSigma
                    varargout{15} = transpose(- 2./(Sigma(ind).^3));
                    if(isempty(varargout{15}))
                        varargout{15} = zeros(1,0);
                    end
                    %dddJ_TdTdRdSigma
                    varargout{16} = transpose(zeros(size(T(ind))));
                    if(isempty(varargout{16}))
                        varargout{16} = zeros(1,0);
                    end
                    %dddJ_TdRdRdSigma
                    varargout{17} = transpose(- 2./(Sigma(ind).^3));
                    if(isempty(varargout{17}))
                        varargout{17} = zeros(1,0);
                    end
                    %dddJ_TdTdSigmadSigma
                    varargout{18} = transpose(6*(T(ind) - Tm(ind))./(Sigma(ind).^4));
                    if(isempty(varargout{18}))
                        varargout{18} = zeros(1,0);
                    end
                    %dddJ_TdRdSigmadSigma
                    varargout{19} = transpose(6*(R(ind))./(Sigma(ind).^4));
                    if(isempty(varargout{19}))
                        varargout{19} = zeros(1,0);
                    end
                    %dddJ_TdSigmadSigmadSigma
                    varargout{20} = transpose(- 12*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^5)) - 12*(((R(ind)).^2)./(Sigma(ind).^5)) + 4./(Sigma(ind).^3));
                    if(isempty(varargout{20}))
                        varargout{20} = zeros(1,0);
                    end
                    if nargout >= 20
                        %ddddJ_TdTdTdTdT
                        varargout{21} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{21}))
                            varargout{21} = zeros(1,0);
                        end
                        %ddddJ_TdTdTdTdR
                        varargout{22} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{22}))
                            varargout{22} = zeros(1,0);
                        end
                        %ddddJ_TdTdTdRdR
                        varargout{23} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{23}))
                            varargout{23} = zeros(1,0);
                        end
                        %ddddJ_TdTdRdRdR
                        varargout{24} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{24}))
                            varargout{24} = zeros(1,0);
                        end
                        %ddddJ_TdRdRdRdR
                        varargout{25} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{25}))
                            varargout{25} = zeros(1,0);
                        end
                        %ddddJ_TdTdTdTdSigma
                        varargout{26} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{26}))
                            varargout{26} = zeros(1,0);
                        end
                        %ddddJ_TdTdTdRdSigma
                        varargout{27} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{27}))
                            varargout{27} = zeros(1,0);
                        end
                        %ddddJ_TdTdRdRdSigma
                        varargout{28} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{28}))
                            varargout{28} = zeros(1,0);
                        end
                        %ddddJ_TdRdRdRdSigma
                        varargout{29} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{29}))
                            varargout{29} = zeros(1,0);
                        end
                        %ddddJ_TdTdTdSigmadSigma
                        varargout{30} = transpose(6./(Sigma(ind).^4));
                        if(isempty(varargout{30}))
                            varargout{30} = zeros(1,0);
                        end
                        %ddddJ_TdTdRdSigmadSigma
                        varargout{31} = transpose(zeros(size(T(ind))));
                        if(isempty(varargout{31}))
                            varargout{31} = zeros(1,0);
                        end
                        %ddddJ_TdRdRdSigmadSigma
                        varargout{32} = transpose(6./(Sigma(ind).^4));
                        if(isempty(varargout{32}))
                            varargout{32} = zeros(1,0);
                        end
                        %ddddJ_TdTdSigmadSigmadSigma
                        varargout{33} = transpose(- 24*((T(ind) - Tm(ind))./(Sigma(ind).^5)));
                        if(isempty(varargout{33}))
                            varargout{33} = zeros(1,0);
                        end
                        %ddddJ_TdRdSigmadSigmadSigma
                        varargout{34} = transpose(- 24*((R(ind))./(Sigma(ind).^5)));
                        if(isempty(varargout{34}))
                            varargout{34} = zeros(1,0);
                        end
                        %ddddJ_TdSigmadSigmadSigmadSigma
                        varargout{35} = transpose(60*(((T(ind) - Tm(ind)).^2)./(Sigma(ind).^6)) + 60*(((R(ind)).^2)./(Sigma(ind).^6)) - 12./(Sigma(ind).^4));
                        if(isempty(varargout{35}))
                            varargout{35} = zeros(1,0);
                        end
                    end
                end
            end
        end
    end
end