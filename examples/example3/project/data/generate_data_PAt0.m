function [Data] = generate_data_PA(Data,Model,s,xi,datafile)
    
    beta = Model.exp{s}.beta(xi);
    [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
    
    Data{s}.PA.beta = beta;
    Data{s}.PA.D = D;
    
    % Simulation of SP approximation
    
    switch(datafile)
        case 'synthetic'
            for i = 1:Model.exp{s}.N
                b_i = sqrtm(D)*randn(size(D,2),1);
                Data{s}.PA.b_i(:,i) = b_i;
                phi_i = Model.exp{s}.phi(beta,b_i);
                [~,~,~,Y] = Model.exp{s}.model(Data{s}.PA.time,phi_i,Data{s}.condition);
                Data{s}.PA.Y(:,:,i) = Y;
                
            end
        otherwise
            num = xlsread(['./project/data/' datafile '.xls']);
            Data{s}.PA.Y(:,1,:) = num(:,2:min(2+Model.exp{s}.N,size(num,2)));
            Data{s}.PA.time = num(:,1);
    end
    
    mdata = mean(Data{s}.PA.Y,3);
    if isfield(Model.exp{s},'PA_post_processing')
        dummy = zeros([size(mdata) size(xi,1)]);
        [mdata,~] = Model.exp{s}.PA_post_processing(mdata,dummy);
    end
    
    sigma = Model.exp{s}.sigma_noise;
    if(size(sigma,1) == size(mdata,1))
        if(size(sigma,2) == 1)
            Sigma = repmat(sigma,[1,size(mdata,2)]);
        elseif(size(sigma,2) == size(mdata,2))
            Sigma = sigma;
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(size(sigma,2) == size(mdata,2))
        if(size(sigma,1) == 1)
            Sigma = repmat(sigma,[size(mdata,1),1]);
        else
            error('Incompatible size of sigma parametrisation!')
        end
    elseif(and(size(sigma,1)==1,size(sigma,2)==1))
        Sigma = repmat(sigma,size(mdata));
    else
        error('Incompatible size of sigma parametrisation!')
    end
    
    Data{s}.PA.m = mdata + Model.exp{s}.noise_on*Sigma.*randn(size(mdata));
    Data{s}.PA.Sigma_m = Sigma;%Model.exp{s}.sigma_noise*ones(size(Data{s}.PA.m));
    
    %     [Data{s}.PA.m,~,~,~] = ...
    %         getSigmaPointApp(@(phi) simulate_RRE_mRNA_transfection_deg_SP(Data{s}.PA.time,phi,Data{s}.condition),...
    %                             Model.exp{s}.A,Model.exp{s}.B,beta,D,dbetadxi,dD_full);
    %     Data{s}.PA.m = Data{s}.PA.m + Model.exp{s}.noise_on*sigma_noise(4)*randn(size(Data{s}.PA.m));
    %     Data{s}.PA.Sigma_m = sigma_noise(4)*ones(size(Data{s}.PA.m));
    
end

