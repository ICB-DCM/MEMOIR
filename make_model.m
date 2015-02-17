function [ Model ] = make_model( Model )
    
    
    if(any([length(Model.param)~=length(Model.fixed_effect),length(Model.param)~=length(Model.random_effect),length(Model.fixed_effect)~=length(Model.random_effect)]))
       error('Size of Model.param, Model.fixed_effect and Model.random_effect do not agree') 
    end
    
    % strip parameters of whitespaces as this will cause problems with generation of symbolic variables
    for j = 1:length(Model.param)
        Model.param{j} = strrep(Model.param{j},' ','_');
    end
    
    
    
    switch(Model.type_D)
        case 'matrix-logarithm'
            n = length(Model.param);
            n_f = sum(Model.fixed_effect);
            n_r = sum(Model.random_effect);
            % initialise symbolic variables
            Model.sym.xi = sym('xi',[n_f+(n_r^2+n_r)/2,1]);
            Model.sym.phi = sym('phi',[sum(arrayfun(@or,Model.fixed_effect,Model.random_effect)),1]);
            Model.sym.beta = sym('beta',[n_f,1]);
            Model.sym.b = sym('b',[n_r,1]);
            Model.sym.delta = sym('delta',[(n_r^2+n_r)/2,1]);
            k = 1;
            k_f = 1;
            k_r = 1;
            
            for j = 1:n;
                if(Model.fixed_effect(j))
                    % construct variables corresponding to fixed effects
                    eval(['syms m_' Model.param{j} ' beta_' Model.param{j} ';']);
                    eval(['Model.sym.xi(k_f) = m_' Model.param{j} ';']);
                    eval(['Model.sym.beta(k_f) = beta_' Model.param{j} ';']);
                    k_f = k_f + 1;
                end
                if(Model.random_effect(j))
                    % construct variables corresponding to radom effects
                    
                    eval(['syms b_' Model.param{j} ';']);
                    ll = find(Model.random_effect(1:k_r+1));
                    for l = 1:length(ll)
                        m = ll(l);
                        eval(['syms C_' Model.param{j} '_' Model.param{m} ' delta_' Model.param{j} '_' Model.param{m} ';']);
                        eval(['Model.sym.xi(sum(Model.fixed_effect)+((k_r-1)^2+(k_r-1))/2+l) = C_' Model.param{j} '_' Model.param{m} ';']);
                        eval(['Model.sym.delta(((k_r-1)^2+(k_r-1))/2+l) = delta_' Model.param{j} '_' Model.param{m} ';'])
                    end
                    eval(['Model.sym.b(k_r) = b_' Model.param{j} ';']);
                    
                    k_r = k_r + 1;
                end
                if(or(Model.fixed_effect(j),Model.random_effect(j)))
                    % construct variables corresponding to mixed effects
                    eval(['syms p_' Model.param{j} ';']);
                    eval(['Model.sym.phi(k) = p_' Model.param{j} ';']);
                    k = k + 1;
                end
            end
            for j = 1:length(Model.random_effect)
                
            end
            
            
        case 'diag-matrix-logarithm'
            n = length(Model.param);
            n_f = sum(Model.fixed_effect);
            n_r = sum(Model.random_effect);
            
            % initialise symbolic variables
            Model.sym.xi = sym('xi',[n_f+n_r,1]);
            Model.sym.phi = sym('phi',[sum(arrayfun(@or,Model.fixed_effect,Model.random_effect)),1]);
            Model.sym.beta = sym('beta',[n_f,1]);
            Model.sym.b = sym('b',[n_r,1]);
            Model.sym.delta = sym('delta',[n_r,1]);
            k = 1;
            k_f = 1;
            k_r = 1;
            
            for j = 1:n;
                if(Model.fixed_effect(j))
                    % construct variables corresponding to fixed effects
                    eval(['syms M_' Model.param{j} ';']);
                    eval(['Model.sym.xi(k_f) = M_' Model.param{j} ';']);
                    eval(['Model.sym.beta(k_f) = M_' Model.param{j} ';']);
                    k_f = k_f + 1;
                end
                if(Model.random_effect(j))
                    % construct variables corresponding to radom effects
                    eval(['syms C_' Model.param{j} ' b_' Model.param{j} ';']);
                    eval(['Model.sym.xi(sum(Model.fixed_effect)+k_r) = C_' Model.param{j} ';']);
                    eval(['Model.sym.b(k_r) = b_' Model.param{j} ';']);
                    eval(['Model.sym.delta(k_r) = C_' Model.param{j} ';']);
                    k_r = k_r + 1;
                end
                if(or(Model.fixed_effect(j),Model.random_effect(j)))
                    % construct variables corresponding to mixed effects
                    eval(['syms p_' Model.param{j} ';']);
                    eval(['Model.sym.phi(k) = p_' Model.param{j} ';']);
                    k = k + 1;
                end
            end
            for j = 1:length(Model.random_effect)
                
            end
            
    end
    
end

