function [Data] = generate_data_SCTLstat(Data,Model,s,xi,datafile)  
    
    beta = Model.exp{s}.beta(xi);
    [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
    
    Data{s}.SCTLstat.beta = beta;
    Data{s}.SCTLstat.D = D;
    
    % Simulation
    switch(datafile)
        case 'synthetic'
            for i = 1:Model.exp{s}.N
                b_i = sqrtm(D)*randn(size(D,2),1);
                Data{s}.SCTLstat.b_i(:,i) = b_i;
                phi_i = Model.exp{s}.phi(beta,b_i);
                options.sensi = 0;
                sol = Model.exp{s}.model(Data{s}.SCTLstat.time,phi_i,Data{s}.condition,options);
                Data{s}.SCTLstat.Y(:,:,i) = sol.y;
            end
        otherwise
            num = xlsread(['./project/data/' datafile '.xls']);
            Data{s}.SCTLstat.Y(:,1,:) = num(:,2:min(2+Model.exp{s}.N,size(num,2)));
            Data{s}.SCTLstat.time = num(:,1);
    end
    
    Z = reshape(Data{s}.SCTLstat.Y,[size(Data{s}.SCTLstat.Y,1)*size(Data{s}.SCTLstat.Y,2),size(Data{s}.SCTLstat.Y,3)]);
    
    mz = mean(Z,2);
    
    sigma = Model.exp{s}.sigma_noise(phi_i);
    if(size(sigma,1) == size(mz,1))
        if(size(sigma,2) == 1)
            Sigma = repmat(sigma,[1,size(mz,2)]);
        elseif(size(sigma,2) == size(mz,2))
            Sigma = sigma;
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == size(mz,2))
        if(size(sigma,1) == 1)
            Sigma = repmat(sigma,[size(mz,1),1]);
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        Sigma = repmat(sigma,size(mz));
    else
        error('Incompatible size of sigma parametrisation!')
    end
    
    Data{s}.SCTLstat.mz = mz + Model.exp{s}.noise_on*Sigma;
    Data{s}.SCTLstat.Sigma_mz = Sigma*ones(size(Z,1),1);
    Data{s}.SCTLstat.Cz = cov(Z') + Model.exp{s}.noise_on*Sigma;
    Data{s}.SCTLstat.Sigma_Cz = Sigma*ones(size(Z,1));

%     % Simulation of SP approximation
%     [~,~,~,Data{s}.SCTLstat.mz,Data{s}.SCTLstat.Cz] = ...
%         getSigmaPointApp(@(phi) simulate_RRE_mRNA_transfection_SP(Data{s}.SCTLstat.time,phi,Data{s}.condition),...
%                             Model.exp{s}.A,Model.exp{s}.B,beta,D,dbetadxi,dD_full);
%     Data{s}.SCTLstat.Sigma_mz = sigma_noise(5)*ones(size(Z,1),1);
%     Data{s}.SCTLstat.Sigma_Cz = sigma_noise(5)*ones(size(Z,1));

end