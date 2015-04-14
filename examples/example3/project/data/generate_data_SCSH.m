function [Data] = generate_data_SCSH(Data,Model,s,xi,datafile)
    
    beta = Model.exp{s}.beta(xi);
    [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
    
    Data{s}.SCSH.beta = beta;
    Data{s}.SCSH.D = D;
    
    % Simulation of SP approximation
    switch(datafile)
        case 'synthetic'
            for i = 1:Model.exp{s}.N
                b_i = sqrtm(D)*randn(size(D,2),1);
                Data{s}.SCSH.b_i(:,i) = b_i;
                phi_i = Model.exp{s}.phi(beta,b_i);
                option_simu.sensi = 0;
                sol = Model.exp{s}.model(Data{s}.SCSH.time,phi_i,Data{s}.condition,option_simu);
                Y = sol.y;
                
                sigma = Model.exp{s}.sigma_noise(phi_i);
                if(size(sigma,1) == size(Y,1))
                    if(size(sigma,2) == 1)
                        Sigma = repmat(sigma,[1,size(Y,2)]);
                    elseif(size(sigma,2) == size(Y,2))
                        Sigma = sigma;
                    else
                        error('Incompatible size of sigma parametrisation!')
                    end
                elseif(size(sigma,2) == size(Y,2))
                    if(size(sigma,1) == 1)
                        Sigma = repmat(sigma,[size(Y,1),1]);
                    else
                        error('Incompatible size of sigma parametrisation!')
                    end
                elseif(and(size(sigma,1)==1,size(sigma,2)==1))
                    Sigma = repmat(sigma,size(Y));
                else
                    error('Incompatible size of sigma parametrisation!')
                end
                
                Data{s}.SCSH.Y(:,:,i) = Y + Model.exp{s}.noise_on*Sigma.*randn(size(Y));
            end
        otherwise
            num = xlsread(['./project/data/' datafile '.xls']);
            Data{s}.SCSH.Y(:,1,:) = num(:,2:min(2+Model.exp{s}.N,size(num,2)));
            Data{s}.SCSH.time = num(:,1);
    end
    
    
    
    %     Data{s}.SCSH.m = mean(Data{s}.SCSH.Y,3) + Model.exp{s}.noise_on*Model.exp{s}.sigma_noise*randn(size(mean(Data{s}.SCSH.Y,3)));
    %     Data{s}.SCSH.Sigma_m = Model.exp{s}.sigma_noise*ones(size(Data{s}.SCSH.m));
    
    for k = 1:length(Data{s}.SCSH.time)
        Data{s}.SCSH.m(k,:) = mean(Data{s}.SCSH.Y(k,:,:),3);
        Data{s}.SCSH.Sigma_m(k,:) = Sigma(k,:);
        Data{s}.SCSH.C(k,:,:) = cov(permute(Data{s}.SCSH.Y(k,:,:),[3,2,1]),1) + Model.exp{s}.noise_on*Sigma(k,:).*randn(size(Data{s}.SCSH.m,2));
        Data{s}.SCSH.Sigma_C(k,:,:) =  sqrt(mean(bsxfun(@minus,Data{s}.SCSH.Y(k,:,:),Data{s}.SCSH.m(k,:)).^4,3)...
                    -(Model.exp{s}.N-3)/(Model.exp{s}.N-1)*diag(Data{s}.SCSH.C(k,:,:)).^2)/sqrt(Model.exp{s}.N);
    end
    
    %Data{s}.SCSH.Sigma_C = Model.exp{s}.sigma_noise*ones(size(Data{s}.SCSH.C));
    
%     [Data{s}.SCSH.m,Data{s}.SCSH.C] = ...
%         getSigmaPointApp(@(phi) simulate_RRE_mRNA_transfection_SP(Data{s}.SCSH.time,phi,Data{s}.condition),...
%         Model.exp{s}.A,Model.exp{s}.B,beta,D,dbetadxi,dD_full);
%     Data{s}.SCSH.m = Data{s}.SCSH.m + Model.exp{s}.noise_on*Model.exp{s}.sigma_noise*randn(size(Data{s}.SCSH.m));
%     Data{s}.SCSH.Sigma_m = Model.exp{s}.sigma_noise*ones(size(Data{s}.SCSH.m));
%     Data{s}.SCSH.C = Data{s}.SCSH.C + Model.exp{s}.noise_on*Model.exp{s}.sigma_noise*randn(size(Data{s}.SCSH.C));
%     Data{s}.SCSH.Sigma_C = Model.exp{s}.sigma_noise*ones(size(Data{s}.SCSH.C));
end