clear
close all
%clc

rng(0);

%% General options

% specify the model name you can add more in the switch definition to
% easily switch between models
model = 'toy';

% specify the datafile here, this will be passed to your specified data_fun
% as argument, this allows you to easily switch between different
% datasets/artificial data
datafile = 'synthetic';

% specify the filename here, this will be used as a name to generate the
% final results file
filename = 'example_1';

switch(model)
    case 'toy'
        % specify the model_fun here this. function should use the true
        % parameters or a guess to generate as basic Model struct aswell as
        % the parameters_MEM struct
        model_fun = @(logbeta,logdelta,logsigma) model_toy(logbeta,logdelta,logsigma);
        % specify the experiment_fun here. this function fills the Model
        % struct with experiments. S is a vector containing the indexes of
        % experiments to consider
        experiment_fun = @(Model,S) experiments_toy(Model,S);
        % specify the data_fun here. this function generates the Data
        % struct by either loading experimental data from the datafile or
        % by generating synthetic data
        data_fun = @(Data,Model,xi,S,datafile) data_toy(Data,Model,xi,S,datafile);
end

% check the date for indexing of results files
date_now = date;
% check the clock for indexing of results files
time = clock;
time_now = [num2str(time(4)) num2str(time(5))];

% generate possible output directories
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
    case 'toy'
        % specify the true/guessed parameters here
        logbeta = log([1;3]);
        logdelta = log([0.7]);
        logsigma = [];
end

xi = [logbeta;logdelta;logsigma];

% Index of experiments
S = [1]; 
% exp_decay
% 1 : Single-cell time lapse

% generate the Model struct
[Model,parameters_MEM] =  model_fun(logbeta,logdelta,logsigma);

% indicate whether integration should be activated for SCTL data
Model.integration = false;
% indicate whether additional penalty terms should be considered for SCTL
% data
Model.penalty = false;
% indicate the model name
Model.name = model;

%% Definition of Experiments

% fill the model struct with experiments
Model = experiment_fun(Model,S);

% compute symbolic derivatives, generate respective m-files and attach them
% to the Model struct
Model = complete_model(Model,S);

%% Data generation

Data = [];

% generate the Data struct
Data = data_fun(Data,Model,xi,S,datafile);

fh_sc = figure;
hold on

%% Mixed Effects Optimization

clc
% Options for multi-start optimization
options.localOptimizerOptions = optimset('algorithm','interior-point',...
    'GradObj','on',...
    'MaxIter',1000,...
    'TolFun',1e-8,...
    'TolX',1e-10,...
    'display','iter',... %        'display','final',...
    'MaxFunEvals',1000*parameters_MEM.number);
options.n_starts = 20;
options.comp_type = 'sequential';
options.mode = 'visual';
options.calc_profiles = true;
options.parameter_index = [3,4];
options = PestoOptions(options);
parameters_MEM.true = xi;

parameters_MEM = getMultiStarts(parameters_MEM,@(theta) logLMEMOIR(theta,Data,Model),options);
parameters_MEM = getProfiles(parameters_MEM,@(theta) logLMEMOIR(theta,Data,Model),options);
save(['./project/results/' model '/' datafile '/' date_now '/MS_MEM/S=[' strrep(num2str(S),'  ','_') ']_' filename '.mat'],...
    'parameters_MEM',...
    'Data',...
    'Model',...
    'options',...
    'S',...
    'xi');