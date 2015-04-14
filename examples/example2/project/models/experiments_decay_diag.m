function Model = experiments_decay_diag(Model,S)

b  = Model.sym.b;
beta = Model.sym.beta;

s = 0;
% Experiment 1
if ismember(1,S)
    s = s + 1;
    Model.exp{s}.N = 10;
    Model.exp{s}.sigma_noise = 1;
    Model.exp{s}.noise_on = 1;
    Model.exp{s}.t = [0:0.01:1]';
    
    Model.exp{s}.ind_phi = [1,2];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2)];
    Model.exp{s}.sym.sigma_noise = sym(0.01);
    Model.exp{s}.sym.sigma_time = sym.empty(0,1);
    
    Model.exp{s}.time_model = 'normal';
    Model.exp{s}.noise_model = 'normal';
    Model.exp{s}.parameter_model = 'normal';
    
    Model.exp{s}.model = @(t,phi,kappa,options) simulate_anal_decay(t,phi,kappa,options);
    
    Model.exp{s}.fh = figure('Name','Decay - SCTL');
    Model.exp{s}.fp = figure('Name','Parameters - SCTL');
    Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCTL(Data,Sim,fh);
end