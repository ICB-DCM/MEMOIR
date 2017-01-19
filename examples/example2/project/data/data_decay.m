function Data = data_decay(Data,Model,xi,S,datafile)
    
    s = 0;
    % Decay experiment: Single-cell time lapse
    if ismember(1,S)
        s = s + 1;
        Data{s}.SCTL.time = Model.exp{s}.t;
        Data{s}.condition = [];
        Data{s}.name = 'Decay - SCTL';
        Data{s}.measurands = {'Signal','Reverse Signal'};
        
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
        Model.exp{s}.plot(Data{s},[],s);
    end

end