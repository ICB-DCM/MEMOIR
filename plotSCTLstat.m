% function fh = plotSCTLstatstat(Data,Sim,fh,options)
function fh = plotSCTLstatstat(varargin)

%% Check and assign inputs
if nargin >= 2
    Data = varargin{1};
    Sim = varargin{2};    
else
    error('Not enough inputs.')
end

% Figure handel
if nargin >= 3
    if ~isempty(varargin{3})
        fh = figure(varargin{3});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Options
options.data.col = 'b';
options.data.ls = '-';
options.data.lw = 1;
options.sim.col = 'r';
options.sim.ls = '--';
options.sim.lw = 1;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% Subplot dimensions
n_y = size(Data.SCTLstat.Y,2);
if ~isempty(Sim)
    nc = 2;
    nr = n_y;
else
    nc = ceil(sqrt(n_y));
    nr = ceil(n_y/nc);
end

%% Visualization: Data and Simulation of individuals
if ~isempty(Sim)
    % Data
    % Means
    subplot(2,2,1); hold off;
    plot(1:length(Data.SCTLstat.mz),Data.SCTLstat.mz,'-o',...
        'linewidth',options.data.lw,...
        'linestyle',options.data.ls,...
        'color',options.data.col); hold on;
    plot(1:length(Data.SCTLstat.mz),Sim.mz,'-x',...
        'linewidth',options.sim.lw,...
        'linestyle',options.sim.ls,...
        'color',options.sim.col); hold on;
    xlabel('(t,y(t))'); ylabel('mean');
    xlim([0.5,length(Data.SCTLstat.mz)+0.5]);

    % Covariances
    subplot(2,2,3); hold off;
    plot(1:length(Data.SCTLstat.Cz(:)),Data.SCTLstat.Cz(:),'-o',...
        'linewidth',options.data.lw,...
        'linestyle',options.data.ls,...
        'color',options.data.col); hold on;
    plot(1:length(Data.SCTLstat.Cz(:)),Sim.Cz(:),'-x',...
        'linewidth',options.sim.lw,...
        'linestyle',options.sim.ls,...
        'color',options.sim.col); hold on;
    xlabel('pairs of (t,y(t))'); ylabel('covariance');
    xlim([0.5,length(Data.SCTLstat.Cz(:))+0.5]);

    % Error
    % Means
    subplot(2,2,2); hold off;
    plot(1:length(Data.SCTLstat.mz),Data.SCTLstat.mz-Sim.mz,'-o',...
        'linewidth',options.error.lw,...
        'linestyle',options.error.ls,...
        'color',options.error.col); hold on;
    xlabel('(t,y(t))'); ylabel('error of mean');
    xlim([0.5,length(Data.SCTLstat.mz)+0.5]);

    % Covariances
    subplot(2,2,4); hold off;
    plot(1:length(Data.SCTLstat.Cz(:)),Data.SCTLstat.Cz(:)-Sim.Cz(:),'-o',...
        'linewidth',options.error.lw,...
        'linestyle',options.error.ls,...
        'color',options.error.col); hold on;
    xlabel('pairs of (t,y(t))'); ylabel('error of covariance');
    xlim([0.5,length(Data.SCTLstat.Cz(:))+0.5]);
end

%% Visualization: Data of individuals
if isempty(Sim)
    % Means
    subplot(2,1,1); hold off;
    plot(1:length(Data.SCTLstat.mz),Data.SCTLstat.mz,'-o',...
        'linewidth',options.data.lw,...
        'linestyle',options.data.ls,...
        'color',options.data.col); hold on;
    xlabel('(t,y(t))'); ylabel('mean');
    xlim([0.5,length(Data.SCTLstat.mz)+0.5]);

    % Covariances
    subplot(2,1,2); hold off;
    plot(1:length(Data.SCTLstat.Cz(:)),Data.SCTLstat.Cz(:),'-o',...
        'linewidth',options.data.lw,...
        'linestyle',options.data.ls,...
        'color',options.data.col); hold on;
    xlabel('pairs of (t,y(t))'); ylabel('covariance');
    xlim([0.5,length(Data.SCTLstat.Cz(:))+0.5]);
end

%%
drawnow

end

