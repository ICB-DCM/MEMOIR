clear
close all
%clc

rng(0);

%% General options

% model = 'transfection_diag';
% model = 'transfection_full';
model = 'transfection_red_diag';
% model = 'transfection_red_full';
% model = 'transfection_logndelay_diag';
% model = 'toy_diag';
% model = 'toy_full';
% model = 'exp_decay';
% datafile = 'Huh7_d2EGFP_130827';
datafile = 'synthetic';

filename = 's=10_woint_hf';

switch(model)
    case 'transfection_diag'
        model_fun = @(logbeta,logdelta,logsigma) model_transfection_diag(logbeta,logdelta,logsigma);
        experiment_fun = @(Model,S) experiments_transfection_diag(Model,S);
        data_fun = @(Data,Model,xi,S,datafile) data_transfection_diag(Data,Model,xi,S,datafile);
    case 'transfection_red_diag'
        model_fun = @(logbeta,logdelta,logsigma) model_transfection_red_diag(logbeta,logdelta,logsigma);
        experiment_fun = @(Model,S) experiments_transfection_red_diag(Model,S);
        data_fun = @(Data,Model,xi,S,datafile) data_transfection_red_diag(Data,Model,xi,S,datafile);
    case 'toy_diag'
        model_fun = @(logbeta,logdelta,logsigma) model_toy_diag(logbeta,logdelta,logsigma);
        experiment_fun = @(Model,S) experiments_toy_diag(Model,S);
        data_fun = @(Data,Model,xi,S,datafile) data_toy_diag(Data,Model,xi,S,datafile);
    case 'toy_full'
        model_fun = @(logbeta,logdelta,logsigma) model_toy_full(logbeta,logdelta,logsigma);
        experiment_fun = @(Model,S) experiments_toy_full(Model,S);
        data_fun = @(Data,Model,xi,S,datafile) data_toy_full(Data,Model,xi,S,datafile);
    case 'exp_decay'
        model_fun = @(logbeta,logdelta,logsigma) model_exp_decay(logbeta,logdelta,logsigma);
        experiment_fun = @(Model,S) experiments_exp_decay(Model,S);
        data_fun = @(Data,Model,xi,S,datafile) data_exp_decay(Data,Model,xi,S,datafile);
end
date_now = date;
time = clock;
time_now = [num2str(time(4)) num2str(time(5))];

comp_SC = 0;
comp_MEM = 1;

eval(['cd ' strrep(mfilename('fullpath'),mfilename('full'),'')])

if(~exist(['./project/results/' model],'dir'))
    mkdir(['./project/results/' model])
    mkdir(['./project/results/' model '/' datafile])
    mkdir(['./project/results/' model '/' datafile '/' date_now '/MS_SC'])
    mkdir(['./project/results/' model '/' datafile '/' date_now '/MS_MEM'])
    mkdir(['./project/results/' model '/' datafile '/' date_now '/PL_MEM'])
else
    if(~exist(['./project/results/' model '/' datafile],'dir'))
        mkdir(['./project/results/' model '/' datafile])
        mkdir(['./project/results/' model '/' datafile '/' date_now '/MS_SC'])
        mkdir(['./project/results/' model '/' datafile '/' date_now '/MS_MEM'])
        mkdir(['./project/results/' model '/' datafile '/' date_now '/PL_MEM'])
    else
        if(~exist(['./project/results/' model '/' datafile '/' date_now '/MS_SC'],'dir'))
            mkdir(['./project/results/' model '/' datafile '/' date_now '/MS_SC'])
        end
        if(~exist(['./project/results/' model '/' datafile '/' date_now '/MS_MEM'],'dir'))
            mkdir(['./project/results/' model '/' datafile '/' date_now '/MS_MEM'])
        end
        if(~exist(['./project/results/' model '/' datafile '/' date_now '/PL_MEM'],'dir'))
            mkdir(['./project/results/' model '/' datafile '/' date_now '/PL_MEM'])
        end
    end
end

%% Definition of mixed-effect model

switch(model)
    case 'transfection_diag'
        % Parameters
        
        logbeta = log([2.5;100;0.5;0.1;25]);
        logdelta = [-4;-4;-2.5;-1.5];
        logsigma = 0;
        % beta = log([2;10;0.2;0.1;0.1]);
        % logdelta = log([0.2;0.2;0.2;0.2]);
    case 'transfection_full'
        % Parameters
        
        logbeta = log([2.5;11;0.5;0.1;25]);
        logdelta = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2];
        % beta = log([2;10;0.2;0.1;0.1]);
        % logdelta = log([0.2;0.2;0.2;0.2]);
    case 'transfection_red_diag'
        % Parameters
        logbeta = log([2.5;100;0.5;0.1]);
        logdelta = [-4;-4;-2.5];
        logsigma = [];
    case 'transfection_red_full'
        % Parameters
        logbeta = log([2.5;100;0.5;0.1]);
        logdelta = [-2;-2;-2;-2;-2;-2];
    case 'transfection_logndelay_diag'
        % Parameters
        logbeta = log([exp(1.5);exp(1e-0);100;0.5;0.1]);
        logdelta = [-4;-3;-4;-2.5;-1.5];
    case 'toy_diag'
        logbeta = log([3;9]);
        logdelta = log([0.1;0.1]);
        logsigma = [];
    case 'toy_full'
        logbeta = log([3;9]);
        logdelta = log([0.1;0.01;0.1]);
        logsigma = [];
    case 'exp_decay'
        logbeta = log([1;3]);
        logdelta = log([0.7]);
        logsigma = [];
end

xi = [logbeta;logdelta;logsigma];

S = [1]; % Index of experiments
% transfection_diag(_full)
% 1 Transfection experiment: Single-cell time lapse
% 2 Degradation experiment: Single-cell time lapse
% 3 Transfection experiment: Population snapshot
% 4 Degradation experiment: Population average
% 5 Transfection experiment: Single-cell time lapse - Statistics
% transfection_red_diag(_full)
% 1 Transfection experiment: Single-cell time lapse
% 2 Transfection experiment: Population snapshot
% 3 Transfection experiment: Single-cell time lapse - Statistics
% toy_diag(_full)
% 1 toy example: Single-cell time lapse

[Model,parameters_MEM] =  model_fun(logbeta,logdelta,logsigma);

Model.integration = false;
Model.penalty = false;
Model.name = model;

%% Definition of Experiments

Model = experiment_fun(Model,S);

% Model.exp{1}.sigma_noise = 10;
% Model.exp{1}.sym.sigma = sym(Model.exp{1}.sigma_noise);

Model = complete_model(Model);

%% Data generation
% Measurement properties

Data = [];

Data = data_fun(Data,Model,xi,S,datafile);



fh_sc = figure;
hold on

if(comp_SC)
    for s = S
        for i = 1:Model.exp{s}.N
            parameters.min = [parameters_MEM.min(Model.exp{s}.ind_beta)];
            parameters.max = [parameters_MEM.max(Model.exp{s}.ind_beta)];
            parameters.number = length(Model.exp{s}.ind_beta);
            parameters.name = [parameters_MEM.name(Model.exp{s}.ind_beta)];
            
            options.fmincon = optimset('algorithm','trust-region-reflective',...
                'GradObj','on',...
                'MaxIter',300,...
                'display','off',...
                'Hessian','user-supplied',...
                'MaxFunEvals',1000*parameters.number);
            options.n_starts = 10;
            options.comp_type = 'sequential';
            options.mode = 'silent';
            options.plot = 'false';
            
            parameters_SC{s,i} = getMultiStarts(parameters,@(theta) logL_SC_w_grad_1(theta,Model,Data,s,i),options);
%             parameters_SC_2{s,i} = getMultiStarts(parameters,@(theta) logL_SC_w_grad_2(theta,Model,Data,s,i),options);
            
            [status,~,~,Y,~,sY] = Model.exp{s}.model(Data{s}.SCTL.time,parameters_SC{s,i}.MS.par(:,1),Data{s}.condition);
            
            figure(fh_sc)
            
            hold on
            
            Ym = Data{s}.SCTL.Y(:,:,i);
            plot(Data{s}.SCTL.time,Ym,'bx','MarkerSize',10,'LineWidth',3)
            plot(Data{s}.SCTL.time,Y,'r--','LineWidth',3)

        end
    end
    save(['./project/results/' model '/' datafile '/' date_now '/MS_SC/S=[' strrep(num2str(S),'  ','_') ']_' filename '.mat'],...
        'parameters_SC',...
        'Data',...
        'Model',...
        'options',...
        'S',...
        'xi');
end



%% Mixed Effects Optimization

if(comp_MEM)
    clc
    % Options for multi-start optimization
    options.fmincon = optimset('algorithm','trust-region-reflective',...
        'GradObj','on',...
        'MaxIter',1000,...
        'TolFun',1e-8,... 
        'TolX',1e-10,...
        'display','iter',... %        'display','final',...
        'MaxFunEvals',1000*parameters_MEM.number,...
        'Hessian','user-supplied');
% options.fmincon = optimset('algorithm','interior-point',...
%     'GradObj','on',...
%     'MaxIter',1000,...
%     'display','iter',...
%     'MaxFunEvals',1000*parameters_MEM.number);
    options.n_starts = 20;
    options.comp_type = 'sequential';
    options.mode = 'visual';
    options.calc_profiles = 'true';
%     parameters_MEM.guess = xi;
    parameters_MEM.true = xi;
    
    % Multi-start optimization
    % logL_CE_w_grad_2(parameters_MEM.guess,Data,Model)
    xi_r = xi+1*rand(size(xi));
%     [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi_r,@(theta) logL_CE_w_grad_2(theta,Data,Model),1e-4,1,2)
%     [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi_r,@(theta) logL_CE_w_grad_2(theta,Data,Model),1e-4,2,3)
    parameters_MEM = getMultiStarts(parameters_MEM,@(theta) logL_CE_w_grad_2(theta,Data,Model),options);
%     parameters_MEM = getProfiles(parameters_MEM,@(theta) logL_CE_w_grad_2(theta,Data,Model),options);
    save(['./project/results/' model '/' datafile '/' date_now '/MS_MEM/S=[' strrep(num2str(S),'  ','_') ']_' filename '.mat'],...
        'parameters_MEM',...
        'Data',...
        'Model',...
        'options',...
        'S',...
        'xi');
%     parameters_PL = getProfiles(parameters_MEM,@(theta) logL_CE_w_grad_2(theta,Data,Model),options);
%     save(['./project/results/' model '/' datafile '/' date_now '/PL_MEM/S=[' strrep(num2str(S),'  ','_') ']_' filename '.mat'],...
%         'parameters_PL',...
%         'Data',...
%         'Model',...
%         'options',...
%         'S',...
%         'xi');
end


