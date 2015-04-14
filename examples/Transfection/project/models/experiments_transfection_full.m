function Model = experiments_transfection_red_full(Model,S)

s = 0;
% Experiment 1
if ismember(1,S)
    s = s + 1;
    Model.exp{s}.N = 10;
    Model.exp{s}.sigma_noise(s) = 1;
    Model.exp{s}.noise_on(s) = 0;
    Model.exp{s}.t = [0:0.2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)+b(4)];
    Model.exp{s}.sym.sigma = sym(1);
    
    Model.exp{s}.model = @(t,phi,kappa) simulate_anal_mRNA_transfection(t,phi,kappa);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCTL');
    Model.exp{s}.fp = figure('Name','Parameters - SCTL');
    Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCTL(Data,Sim,fh);
end

% Experiment 2
if ismember(2,S)
    s = s + 1;
    Model.exp{s}.N = 1000;
    Model.exp{s}.sigma_noise(s) = 1;
    Model.exp{s}.noise_on(s) = 0;
    Model.exp{s}.t = [0:0.2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)+b(4)];
    Model.exp{s}.sym.sigma = sym(1);
    
    Model.exp{s}.model = @(t,phi,kappa) simulate_anal_mRNA_transfection(t,phi,kappa);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCSH');
    Model.exp{s}.fp = figure('Name','Parameters - SCTL');
    Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCSH(Data,Sim,fh);
end

% Experiment 3
if ismember(3,S)
    s = s + 1;
    Model.exp{s}.N = 1000;
    Model.exp{s}.sigma_noise(s) = 1;
    Model.exp{s}.noise_on(s) = 0;
    Model.exp{s}.t = [0:2:24]';
    
    Model.exp{s}.ind_phi = [1,2,3,4];
    Model.exp{s}.sym.phi = [beta(1)+b(1);
        beta(2)+b(2);
        beta(3)+b(3);
        beta(4)+b(4)];
    Model.exp{s}.sym.sigma = sym(1);
    
    Model.exp{s}.model = @(t,phi,kappa) simulate_anal_mRNA_transfection(t,phi,kappa);
    
    Model.exp{s}.fh = figure('Name','Transfection - SCTL');
    Model.exp{s}.fp = figure('Name','Parameters - SCTL');
    Model.exp{s}.fl = figure('Name','Likelihood Contribution - SCTL');
    Model.exp{s}.plot = @(Data,Sim,fh) plotSCTLstat(Data,Sim,fh);
end
