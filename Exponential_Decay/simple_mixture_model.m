clear all;
close all;
clc;

%% OPTIONS
alpha = 0.99; % confidence level

%% Data generation
% # individuals
n = 20; 

% true parameters
x0_true = 1;
mu_d_true = 1;
sigma_d_true = 0.7;
% sigma_noise_true = 0.1; % low signal-to-noise ratio
% sigma_noise_true = 0.3; % high signal-to-noise ratio
sigma_noise_true = 20; % very very high signal-to-noise ratio
d_true = random('logn',mu_d_true,sigma_d_true,[n,1]);
theta_true = [log10(x0_true);mu_d_true;log10(sigma_d_true);log10(sigma_noise_true);log10(d_true)];

% measurements
t = [0.2:0.8:1]';
D = zeros(length(t),n);
for i = 1:n
    D(:,i) = x0_true*exp(-d_true(i)*t);
end
D = D.*random('logn',0,sigma_noise_true,[length(t),n]);

% plot
figure('name','Measurement data');
leg_str = {};
C = colormap(jet(n));
for i = 1:n
    leg_h(i) = plot(t,D(:,i),'*-','color',C(i,:)); hold on;
    leg_str{end+1} = ['cell ' num2str(i,'%d')];
end
legend(leg_h,leg_str);
xlabel('time');
ylabel('observable');
set(gca,'yscale','log');

% %% Estimation - individual
% % Parameters
% parameters.min = -2*ones(4+n,1);
% parameters.max =  2*ones(4+n,1);
% parameters.name = {'log_{10}(x_0)';'\mu_d';'log_{10}(\sigma_d)';'log_{10}(\sigma_{noise})'};
% for i = 1:n
%     parameters.name{end+1} = ['log_{10}(d_{' num2str(i,'%d') '})'];
% end
% parameters.number = length(parameters.min);
% 
% % Options for multi-start optimization
% options.fmincon = optimset('algorithm','interior-point',...
%                            'GradObj','on',...
%                            'MaxIter',3000,...
%                            'MaxFunEvals',3000*5);
% options.n_starts = 20;
% options.plot = 'true';
% options.fh = figure('name','Result of multi-start optimization');
% 
% % Log-likelihood function
% switch options.fmincon.GradObj
%     case 'off'
%         logL = @(theta,options) logL_MEM_wo_grad(theta,t,D,options);
%     case 'on'
%         logL = @(theta,options) logL_MEM_w_grad(theta,t,D,options);
% end
% 
% % Multi-start optimization
% parameters = optimizeMultiStart(parameters,logL,options);
% theta_MAP = parameters.MS.MAP.par;
% x0_MAP = 10.^parameters.MS.MAP.par(1);
% mu_d_MAP = parameters.MS.MAP.par(2);
% sigma_d_MAP = 10.^parameters.MS.MAP.par(3);
% sigma_noise_MAP = 10.^parameters.MS.MAP.par(4);
% d_MAP = 10.^parameters.MS.MAP.par(5:end);
% 
% % % Profile likelihood calculation
% % options_PL.fmincon = options.fmincon;
% % options_PL.fh = figure('name','Result of profile likelihood calculation');
% % options_PL.parameter_index = 1:4; % profile calculation only for population parameters
% % parameters = computeProfiles(parameters,logL,options_PL);
% 
% % Confidence intervals
% parameters = getConfidenceIntervals(parameters,alpha);
% 
% % Model comparison
% parameters.AIC = 2*parameters.number - 2*parameters.MS.MAP.logPost;
% parameters.BIC = parameters.number*log(n*length(t)) - 2*parameters.MS.MAP.logPost;
% 
% %% Visulalization - individual
% figure('name','Fit of individual datasets');
% % Initalization
% lh = [];
% t_sim = linspace(0,max(t),30)';
% % Number of rows and columns of figure
% num_cols = ceil(sqrt(n));
% num_rows = ceil(n/num_cols);
% % Loop over individuals
% for i = 1:n
%     subplot(num_rows,num_cols,i);
%     % Simulation
%     y = x0_MAP*exp(-d_MAP(i)*t_sim);
%     y_min = y*exp(-3*sigma_noise_MAP);
%     y_max = y*exp(+3*sigma_noise_MAP);
%     % Plot
%     lh(3) = fill(t_sim([1:end,end:-1:1]),[y_min;y_max(end:-1:1)],0.7*[1,1,1]); hold on;
%     plot(t_sim,y_min,'-','color',0.4*[1,1,1],'linewidth',1);
%     plot(t_sim,y_max,'-','color',0.4*[1,1,1],'linewidth',1);
%     lh(2) = plot(t_sim,y,'-','color',0*[1,1,1],'linewidth',1.5);
%     lh(1) = plot(t,D(:,i),'r*','linewidth',2);
%     % Legend
%     if i == 1
%         legend(lh,{'data','model','model +/- 3 \sigma'},'Location','NorthEast');
%     end
%     xlabel('time');
%     ylabel('observable');
%     title(['cell ' num2str(i,'%d')]);
%     set(gca,'yscale','log');
% end

%% Estimation - marginalized
% Parameters
parameters_m.min = -2*ones(4,1);
parameters_m.max =  2*ones(4,1);
parameters_m.name = {'log_{10}(x_0)';'\mu_d';'log_{10}(\sigma_d)';'log_{10}(\sigma_{noise})'};
parameters_m.number = length(parameters_m.min);

% Options for multi-start optimization
options_m.fmincon = optimset('algorithm','interior-point',...
                             'GradObj','on',...
                             'MaxIter',3000,...
                             'MaxFunEvals',3000*5);
options_m.n_starts = 20;
options_m.plot = 'true';
options_m.fh = figure('name','Result of multi-start otpimization: marginal');

% Log-likelihood function
switch options_m.fmincon.GradObj
    case 'off'
        logL = @(theta,options) logL_MEM_m_wo_grad(theta,t,D,options);
    case 'on'
        logL = @(theta,options) logL_MEM_m_w_grad(theta,t,D,options);
end

% Multi-start optimization
parameters_m = optimizeMultiStart(parameters_m,logL,options_m);
theta_m_MAP = parameters_m.MS.MAP.par;
x0_m_MAP = 10.^parameters_m.MS.MAP.par(1);
mu_d_m_MAP = parameters_m.MS.MAP.par(2);
sigma_d_m_MAP = 10.^parameters_m.MS.MAP.par(3);
sigma_noise_m_MAP = 10.^parameters_m.MS.MAP.par(4);

% % Profile likelihood calculation
% options_m_PL.fmincon = options_m.fmincon;
% options_m_PL.fh = figure('name','Result of profile likelihood calculation: marginal');
% options_m_PL.parameter_index = 1:2; % profile calculation only for population parameters
% parameters_m = computeProfiles(parameters_m,logL,options_m_PL);

% Confidence intervals
% parameters_m = getConfidenceIntervals(parameters_m,alpha);

% Model comparison
parameters_m.AIC = 2*parameters_m.number - 2*parameters_m.MS.MAP.logPost;
parameters_m.BIC = parameters_m.number*log(n*length(t)) - 2*parameters_m.MS.MAP.logPost;

%% Estimation - marginalized: individual parameters
% Parameters
parameters_m_ind.min = -2;
parameters_m_ind.max =  2;
parameters_m_ind.name = {'log_{10}(d)'};
parameters_m_ind.number = 1;
parameters_m_ind.guess = log10(exp(mu_d_m_MAP));

% Options for multi-start optimization
options_m_ind.fmincon = optimset('algorithm','interior-point',...
                                 'GradObj','on',...
                                 'MaxIter',3000,...
                                 'MaxFunEvals',3000*5);
options_m_ind.n_starts = 1;
options_m_ind.plot = 'false';% 'true';
options_m_ind.mode = 'silent';% 'normal';
% options_m_ind.plot = 'true';
% options_m_ind.mode = 'normal';
options_m_ind.logPost_options.grad_ind = 5;

% Loop: single individuals
log10_d_m_MAP = zeros(n,1);
d_m_MAP = zeros(n,1);
for i = 1:n
    % Log-likelihood function
    switch options_m_ind.fmincon.GradObj
        case 'off'
            logL = @(theta,options) logL_MEM_wo_grad([log10(x0_m_MAP),...
                                                      mu_d_m_MAP,...
                                                      log10(sigma_d_m_MAP),...
                                                      log10(sigma_noise_m_MAP),...
                                                      theta],t,D(:,i),options);
        case 'on'
            logL = @(theta,options) logL_MEM_w_grad([log10(x0_m_MAP),...
                                                     mu_d_m_MAP,...
                                                     log10(sigma_d_m_MAP),...
                                                     log10(sigma_noise_m_MAP),...
                                                     theta],t,D(:,i),options);
    end

    % Multi-start optimization
    [parameters_m_ind,fh_MS] = optimizeMultiStart(parameters_m_ind,logL,options_m_ind);

    % Assignment
    log10_d_m_MAP(i) = parameters_m_ind.MS.MAP.par;
    d_m_MAP(i) = 10.^parameters_m_ind.MS.MAP.par;
end

% Construct full parameter vector
theta_m_MAP = [theta_m_MAP;log10_d_m_MAP];

%% Visulalization - marginal
figure('name','Fit of individual datasets: marginal');
% Initalization
lh = [];
t_sim = linspace(0,max(t),30)';
% Number of rows and columns of figure
num_cols = ceil(sqrt(n));
num_rows = ceil(n/num_cols);
% Loop over individuals
for i = 1:n
    subplot(num_rows,num_cols,i);
    % Simulation
    y = x0_m_MAP*exp(-d_m_MAP(i)*t_sim);
    y_min = y*exp(-3*sigma_noise_m_MAP);
    y_max = y*exp(+3*sigma_noise_m_MAP);
    % Plot
    lh(3) = fill(t_sim([1:end,end:-1:1]),[y_min;y_max(end:-1:1)],0.7*[1,1,1]); hold on;
    plot(t_sim,y_min,'-','color',0.4*[1,1,1],'linewidth',1);
    plot(t_sim,y_max,'-','color',0.4*[1,1,1],'linewidth',1);
    lh(2) = plot(t_sim,y,'-','color',0*[1,1,1],'linewidth',1.5);
    lh(1) = plot(t,D(:,i),'r*','linewidth',2);
    % Legend
    if i == 1
        legend(lh,{'data','model','model +/- 3 \sigma'},'Location','NorthEast');
    end
    xlabel('time');
    ylabel('observable');
    title(['cell ' num2str(i,'%d')]);
    set(gca,'yscale','log');
end

% %% Estimation - mean
% % Parameters
% parameters_mean.min = -2*ones(3,1);
% parameters_mean.max =  2*ones(3,1);
% parameters_mean.name = {'log_{10}(x_0)';'log_{10}(d)';'log_{10}(\sigma_{noise})'};
% parameters_mean.number = 3;
% 
% % Options for multi-start optimization
% options_mean.fmincon = optimset('algorithm','interior-point',...
%                                 'GradObj','on',...
%                                 'MaxIter',3000,...
%                                 'MaxFunEvals',3000*5);
% options_mean.n_starts = 20;
% options_mean.plot = 'true';
% options_mean.fh = figure('name','Result of multi-start optimization: mean');
% 
% % Log-likelihood function
% switch options_mean.fmincon.GradObj
%     case 'off'
%         logL = @(theta,options) logL_wo_grad(theta,t,D,options);
%     case 'on'
%         logL = @(theta,options) logL_w_grad(theta,t,D,options);
% end
% 
% % Multi-start optimization
% parameters_mean = optimizeMultiStart(parameters_mean,logL,options_mean);
% theta_mean_MAP = parameters_mean.MS.MAP.par;
% x0_mean_MAP = 10.^parameters_mean.MS.MAP.par(1);
% d_mean_MAP = 10.^parameters_mean.MS.MAP.par(2);
% sigma_noise_mean_MAP = 10.^parameters_mean.MS.MAP.par(3);
% 
% % % Profile likelihood calculation
% % options_mean_PL.fmincon = options_mean.fmincon;
% % options_mean_PL.fh = figure('name','Result of profile likelihood calculation: mean');
% % options_mean_PL.parameter_index = 1:3; % profile calculation only for population parameters
% % parameters_mean = computeProfiles(parameters_mean,logL,options_mean_PL);
% 
% % Confidence intervals
% %parameters_mean = getConfidenceIntervals(parameters_mean,alpha);
% 
% % Model comparison
% parameters_mean.AIC = 2*parameters_mean.number - 2*parameters_mean.MS.MAP.logPost;
% parameters_mean.BIC = parameters_mean.number*log(n*length(t)) - 2*parameters_mean.MS.MAP.logPost;
% 
% %% Visulalization - mean
% figure('name','Fit of individual datasets: mean');
% % Initalization
% lh = [];
% t_sim = linspace(0,max(t),30)';
% % Number of rows and columns of figure
% num_cols = ceil(sqrt(n));
% num_rows = ceil(n/num_cols);
% % Loop over individuals
% for i = 1:n
%     subplot(num_rows,num_cols,i);
%     % Simulation
%     y = x0_MAP*exp(-d_mean_MAP*t_sim);
%     y_min = y*exp(-3*sigma_noise_mean_MAP);
%     y_max = y*exp(+3*sigma_noise_mean_MAP);
%     % Plot
%     lh(3) = fill(t_sim([1:end,end:-1:1]),[y_min;y_max(end:-1:1)],0.7*[1,1,1]); hold on;
%     plot(t_sim,y_min,'-','color',0.4*[1,1,1],'linewidth',1);
%     plot(t_sim,y_max,'-','color',0.4*[1,1,1],'linewidth',1);
%     lh(2) = plot(t_sim,y,'-','color',0*[1,1,1],'linewidth',1.5);
%     lh(1) = plot(t,D(:,i),'r*','linewidth',2);
%     % Legend
%     if i == 1
%         legend(lh,{'data','model','model +/- 3 \sigma'},'Location','NorthEast');
%     end
%     xlabel('time');
%     ylabel('observable');
%     title(['cell ' num2str(i,'%d')]);
%     set(gca,'yscale','log');
% end
% 
% %% 
% [theta_true,theta_MAP,theta_m_MAP]
% 
% [parameters.MS.MAP.logPost,parameters_m.MS.MAP.logPost,parameters_mean.MS.MAP.logPost;
%  parameters.AIC,parameters_m.AIC,parameters_mean.AIC;
%  parameters.BIC,parameters_m.BIC,parameters_mean.BIC]
% 
% save all
