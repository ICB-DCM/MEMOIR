function Data = data_transfection_red_diag(Data,Model,xi,S,datafile)
    
    s = 0;
    % Transfection experiment: Single-cell time lapse
    if ismember(1,S)
        s = s + 1;
        Data{s}.SCTL.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCTL';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
        n_b = length(Model.exp{s}.ind_b);
        
        if(strcmp(datafile,'synthetic'))
            Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
            
            for i = 1:Model.exp{s}.N
                Data{s}.b_i(:,i) = sqrtm(D)*randn(size(D,2),1);
            end
        end
        
        Data = generate_data_SCTL(Data,Model,s,xi,datafile);
        
        % if no sigma is defined, define data-dependant sigma
        if(not(isfield(Model.exp{s},'sigma')))
            Model.exp{s}.sigma = @(xi)  Data{s}.SCTL.Sigma_Y;
            Model.exp{s}.dsigmadxi = @(xi)  zeros(size(Data{s}.SCTL.Sigma_Y));
        end
        
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
        
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
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
        
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
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
    
    % Transfection experiment: Single-cell time lapse
    if ismember(4,S)
        s = s + 1;
        Data{s}.SCTL.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCTL';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
        n_b = length(Model.exp{s}.ind_b);
        
        if(strcmp(datafile,'synthetic'))
            Data{s}.b_i = NaN(n_b,Model.exp{s}.N);
            
            for i = 1:Model.exp{s}.N
                Data{s}.b_i(:,i) = sqrtm(D)*randn(size(D,2),1);
            end
        end
        
        Data = generate_data_SCTL(Data,Model,s,xi,datafile);
        
        % if no sigma is defined, define data-dependant sigma
        if(not(isfield(Model.exp{s},'sigma')))
            Model.exp{s}.sigma = @(xi)  Data{s}.SCTL.Sigma_Y;
            Model.exp{s}.dsigmadxi = @(xi)  zeros(size(Data{s}.SCTL.Sigma_Y));
        end
        
        % Visualization
        Model.exp{s}.plot(Data{s},[],Model.exp{s}.fh);
    end
    
    % Transfection experiment: Population snapshot
    if ismember(5,S)
        s = s + 1;
        Data{s}.SCSH.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCSH';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
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
    if ismember(6,S)
        s = s + 1;
        Data{s}.SCTLstat.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Transfection - SCTLstat';
        Data{s}.measurands = {'GFP'};
        
        [D,~,~] = xi2D(Model.exp{s}.delta(xi),Model.type_D);
        
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