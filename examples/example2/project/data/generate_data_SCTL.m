function [Data] = generate_data_SCTL(Data,Model,s,xi,datafile)
    
    beta = Model.exp{s}.beta(xi);
    [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
    
    Data{s}.SCTL.beta = beta;
    Data{s}.SCTL.D = D;
    
    % Simulation
    
    switch(datafile)
        case 'synthetic'
            for i = 1:Model.exp{s}.N
                b_i = sqrtm(D)*randn(size(D,2),1);
                Data{s}.SCTL.b_i(:,i) = b_i;
                phi_i = Model.exp{s}.phi(beta,b_i);
                option_simu.sensi = 0;
                sol = Model.exp{s}.model(Data{s}.SCTL.time,phi_i,Data{s}.condition,option_simu);
                Y = sol.y;
                if(size(Model.exp{s}.sigma_noise(phi_i),1)==size(Y,1))
                    Data{s}.SCTL.Y(:,:,i) = Y + Model.exp{s}.noise_on*repmat(Model.exp{s}.sigma_noise(phi_i),[length(Data{s}.SCTL.time),1]).*randn(size(Y));
                    Data{s}.SCTL.Sigma_Y(:,:,i) = repmat(Model.exp{s}.sigma_noise(phi_i),[length(Data{s}.SCTL.time),1]);
                else
                    Data{s}.SCTL.Y(:,:,i) = Y + Model.exp{s}.noise_on*Model.exp{s}.sigma_noise(phi_i)*randn(size(Y));
                    Data{s}.SCTL.Sigma_Y(:,:,i) = Model.exp{s}.sigma_noise(phi_i)*ones(size(Y));
                end
            end
        otherwise
            num = xlsread(['./project/data/' datafile '.xls']);
            Data{s}.SCTL.Y(:,1,:) = num(:,2:min(2+Model.exp{s}.N,size(num,2)));
            Data{s}.SCTL.time = num(:,1);
            Data{s}.SCTL.Sigma_Y = ones(size(Data{s}.SCTL.Y));
    end
    
    
end