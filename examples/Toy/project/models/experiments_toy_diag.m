function Model = experiments_toy_diag(Model,S)

    b  = Model.sym.b;
    beta = Model.sym.beta;

s = 0;
% Experiment 1
if ismember(1,S)
    s = s + 1;
    Model.exp{s}.N = 10;
    Model.exp{s}.sigma_noise = 0.1;
    Model.exp{s}.noise_on = 0;
    Model.exp{s}.t = [0:0.01:1]';
    
    Model.exp{s}.ind_phi = [1,2];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2)];
    Model.exp{s}.sym.sigma = sym(Model.exp{s}.sigma_noise);
    
    Model.exp{s}.noise_model = 'normal';
    Model.exp{s}.parameter_model = 'normal';
    
    Model.exp{s}.model = @(t,phi,kappa) nonlinfun(t,phi,kappa);
    
    Model.exp{s}.fh = figure('Name','Toy - SCTL');
%     Model.exp{s}.fp = figure('Name','Parameters - SCTL');
%     Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCTL(Data,Sim,fh);
end