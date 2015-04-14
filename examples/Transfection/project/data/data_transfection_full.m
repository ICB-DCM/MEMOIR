function Data = data_transfection_red_full(Data,Model,xi,S,datafile)
    
    s = 0;
    % Transfection experiment: Single-cell time lapse
    if ismember(1,S)
        s = s + 1;
        Data{s}.SCTL.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCTL';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(xi(Model.exp{s}.ind_D),Model.type_D);
        
        n_b = length(Model.exp{s}.ind_b);
        
        if(strcmp(datafile,'synthetic'))
            Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
            
            for i = 1:Model.exp{s}.N
                Data{s}.b_i(:,i) = sqrtm(D)*randn(size(D,2),1);
            end
        end
        
        Data = generate_data_SCTL(Data,Model,s,xi,datafile);
        
        % Visualization
        Model.exp{s}.plot(Data{s},[],Model.exp{s}.fh);
    end
    
    % Transfection experiment: Population snapshot
    if ismember(2,S)
        s = s + 1;
        Data{s}.SCSH.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCSH';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(xi(Model.exp{s}.ind_D),Model.type_D);
        
        n_b = length(Model.exp{s}.ind_b);
        
        if(strcmp(datafile,'synthetic'))
            Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
            
            for i = 1:Model.exp{s}.N
                Data{s}.b_i(:,i) = sqrtm(D)*randn(size(D,2),1);
            end
        end
        
        Data = generate_data_SCSH(Data,Model,s,xi,datafile);
        
        % Visualization
        Model.exp{s}.plot(Data{s},[],Model.exp{s}.fh);
    end
    
    % Transfection experiment: Single-cell time lapse - Statistics
    if ismember(3,S)
        s = s + 1;
        Data{s}.SCTLstat.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCTLstat';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(xi(Model.exp{s}.ind_D),Model.type_D);
        
        n_b = length(Model.exp{s}.ind_b);
        
        Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
        
        if(strcmp(datafile,'synthetic'))
            Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
            
            for i = 1:Model.exp{s}.N
                Data{s}.b_i(:,i) = sqrtm(D)*randn(size(D,2),1);
            end
        end
        
        Data = generate_data_SCTLstat(Data,Model,s,xi,datafile);
        
        % Visualization
        Model.exp{s}.plot(Data{s},[],Model.exp{s}.fh);
    end
    
end