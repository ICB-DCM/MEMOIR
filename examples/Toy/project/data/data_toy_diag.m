function Data = data_toy_diag(Data,Model,xi,S,datafile)
    
    s = 0;
    % toy example: Single-cell time lapse
    if ismember(1,S)
        s = s + 1;
        Data{s}.SCTL.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'toy example - SCTL';
        Data{s}.measurands = {'x1','x2'};
        
        beta = xi(Model.exp{s}.ind_beta);
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
        n_b = length(Model.exp{s}.ind_b);
        
        Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
        
        for i = 1:Model.exp{s}.N
            Data{s}.b_i(:,i) = sqrtm(D)*randn(size(D,2),1);
        end
        
        Data = generate_data_SCTL(Data,Model,s,xi,datafile);
        
        % Visualization
        Model.exp{s}.plot(Data{s},[],Model.exp{s}.fh);
    end
    
end