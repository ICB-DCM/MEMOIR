function Model = experiments_transfection(Model,S)
    
    b  = Model.sym.b;
    beta = Model.sym.beta;
    
s = 0;
% Experiment 1
if ismember(1,S)
    s = s + 1;
    Model.exp{s}.N = 20;
    Model.exp{s}.sigma_noise = 10;
    Model.exp{s}.sigma_time = 1;
    Model.exp{s}.noise_on = 1;
    Model.exp{s}.t = [0:0.2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)];
    Model.exp{s}.sym.sigma_noise = sym(10);
    Model.exp{s}.sym.sigma_time = sym(1);
    
    Model.exp{s}.noise_model = 'normal';
    Model.exp{s}.parameter_model = 'normal';
    Model.exp{s}.time_model = 'normal';
    
    Model.exp{s}.model = @(t,phi,kappa,options) simulate_anal_mRNA_transfection(t,phi,kappa,options);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCTL');
    Model.exp{s}.fp = figure('Name','Parameters - SCTL');
    Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCTL(Data,Sim,fh);
end

% Experiment 2
if ismember(2,S)
    s = s + 1;
    Model.exp{s}.N = 100;
    Model.exp{s}.sigma_noise = 10;
    Model.exp{s}.sigma_time = 1;
    Model.exp{s}.noise_on = 1;
    Model.exp{s}.t = [0:0.2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)];
    Model.exp{s}.sym.sigma_noise = sym(10);
    Model.exp{s}.sym.sigma_time = sym(1);
    
    Model.exp{s}.noise_model = 'normal';
    Model.exp{s}.parameter_model = 'normal';
    
    Model.exp{s}.model = @(t,phi,kappa,options) simulate_anal_mRNA_transfection(t,phi,kappa,options);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCSH');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCSH(Data,Sim,fh);
end

if ismember(3,S)
    s = s + 1;
    Model.exp{s}.N = 20;
    Model.exp{s}.sigma_noise = 20;
    Model.exp{s}.sigma_time = 1;
    Model.exp{s}.noise_on = 1;
    Model.exp{s}.t = [0:0.2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)];
    Model.exp{s}.sym.sigma_noise = sym(20);
    Model.exp{s}.sym.sigma_time = sym(1);
    
    Model.exp{s}.noise_model = 'normal';
    Model.exp{s}.parameter_model = 'normal';
    Model.exp{s}.time_model = 'normal';
    
    Model.exp{s}.model = @(t,phi,kappa,options) simulate_anal_mRNA_transfection(t,phi,kappa,options);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCTL');
    Model.exp{s}.fp = figure('Name','Parameters - SCTL');
    Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCTL(Data,Sim,fh);
end

% Experiment 2
if ismember(5,S)
    s = s + 1;
    Model.exp{s}.N = 100;
    Model.exp{s}.sigma_noise = 20;
    Model.exp{s}.sigma_time = 1;
    Model.exp{s}.noise_on = 1;
    Model.exp{s}.t = [0:0.2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)];
    Model.exp{s}.sym.sigma_noise = sym(20);
    Model.exp{s}.sym.sigma_time = sym(1);
    
    Model.exp{s}.noise_model = 'normal';
    Model.exp{s}.parameter_model = 'normal';
    
    Model.exp{s}.model = @(t,phi,kappa,options) simulate_anal_mRNA_transfection(t,phi,kappa,options);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCSH');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCSH(Data,Sim,fh);
end