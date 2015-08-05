% function plotSCTLstatstat(Data,Sim,fh,options)
function plotSCTLstat(varargin)

persistent fh

%% Check and assign inputs
if nargin >= 2
    Data = varargin{1};
    Sim = varargin{2};    
else
    error('Not enough inputs.')
end

% Figure handle
s = varargin{3};
if(isempty(fh))
    fh(s) = figure;
else
    if length(fh) < s
        fh(s) = figure;
    elseif(isempty(fh(s)))
        fh(s) = figure;
    end
end
figure(fh(s));

% Options
options.data.col = 'b';
options.data.ls = '-';
options.data.lw = 2;
options.sim.col = 'r';
options.sim.ls = '-';
options.sim.lw = 2;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

% %% Subplot dimensions
% n_y = size(Data.SCTLstat.Y,2);
% if ~isempty(Sim)
%     nc = 2;
%     nr = n_y;
% else
%     nc = ceil(sqrt(n_y));
%     nr = ceil(n_y/nc);
% end

% %% Visualization: Data and Simulation
% if ~isempty(Sim)
%     % Data
%     % Means
%     subplot(2,2,1); hold off;
%     plot(1:length(Data.SCTLstat.mz),Data.SCTLstat.mz,'-o',...
%         'linewidth',options.data.lw,...
%         'linestyle',options.data.ls,...
%         'color',options.data.col); hold on;
%     plot(1:length(Data.SCTLstat.mz),Sim.mz,'-x',...
%         'linewidth',options.sim.lw,...
%         'linestyle',options.sim.ls,...
%         'color',options.sim.col); hold on;
%     xlabel('(t,y(t))'); ylabel('mean');
%     xlim([0.5,length(Data.SCTLstat.mz)+0.5]);
% 
%     % Covariances
%     subplot(2,2,3); hold off;
%     plot(1:length(Data.SCTLstat.Cz(:)),Data.SCTLstat.Cz(:),'-o',...
%         'linewidth',options.data.lw,...
%         'linestyle',options.data.ls,...
%         'color',options.data.col); hold on;
%     plot(1:length(Data.SCTLstat.Cz(:)),Sim.Cz(:),'-x',...
%         'linewidth',options.sim.lw,...
%         'linestyle',options.sim.ls,...
%         'color',options.sim.col); hold on;
%     xlabel('pairs of (t,y(t))'); ylabel('covariance');
%     xlim([0.5,length(Data.SCTLstat.Cz(:))+0.5]);
% 
%     % Error
%     % Means
%     subplot(2,2,2); hold off;
%     plot(1:length(Data.SCTLstat.mz),Data.SCTLstat.mz-Sim.mz,'-o',...
%         'linewidth',options.error.lw,...
%         'linestyle',options.error.ls,...
%         'color',options.error.col); hold on;
%     xlabel('(t,y(t))'); ylabel('error of mean');
%     xlim([0.5,length(Data.SCTLstat.mz)+0.5]);
% 
%     % Covariances
%     subplot(2,2,4); hold off;
%     plot(1:length(Data.SCTLstat.Cz(:)),Data.SCTLstat.Cz(:)-Sim.Cz(:),'-o',...
%         'linewidth',options.error.lw,...
%         'linestyle',options.error.ls,...
%         'color',options.error.col); hold on;
%     xlabel('pairs of (t,y(t))'); ylabel('error of covariance');
%     xlim([0.5,length(Data.SCTLstat.Cz(:))+0.5]);
% end
% 
% %% Visualization: Data
% if isempty(Sim)
%     % Means
%     subplot(2,1,1); hold off;
%     plot(1:length(Data.SCTLstat.mz),Data.SCTLstat.mz,'-o',...
%         'linewidth',options.data.lw,...
%         'linestyle',options.data.ls,...
%         'color',options.data.col); hold on;
%     xlabel('(t,y(t))'); ylabel('mean');
%     xlim([0.5,length(Data.SCTLstat.mz)+0.5]);
% 
%     % Covariances
%     subplot(2,1,2); hold off;
%     plot(1:length(Data.SCTLstat.Cz(:)),Data.SCTLstat.Cz(:),'-o',...
%         'linewidth',options.data.lw,...
%         'linestyle',options.data.ls,...
%         'color',options.data.col); hold on;
%     xlabel('pairs of (t,y(t))'); ylabel('covariance');
%     xlim([0.5,length(Data.SCTLstat.Cz(:))+0.5]);
% end

%% Visualization: Data and Simulation
if ~isempty(Sim)
    % Number of outputs
    ny = size(Data.SCTLstat.mz,2);
    
    % Means
    for i = 1:ny
        subplot(ny+1,2*ny,[2*i-1,2*i]); hold off;
        plot(Data.SCTLstat.time,Data.SCTLstat.mz(:,i),'-',...
            'linewidth',options.data.lw,...
            'linestyle',options.data.ls,...
            'color',options.data.col); hold on;
        plot(Data.SCTLstat.time,Data.SCTLstat.mz(:,i)-3*Data.SCTLstat.Sigma_mz(:,i),'-',...
            'linewidth',options.error.lw,...
            'linestyle',options.error.ls,...
            'color',options.error.col); hold on;
        plot(Data.SCTLstat.time,Data.SCTLstat.mz(:,i)+3*Data.SCTLstat.Sigma_mz(:,i),'-',...
            'linewidth',options.error.lw,...
            'linestyle',options.error.ls,...
            'color',options.error.col); hold on;
        plot(Data.SCTLstat.time,Sim.mz(:,i),'-',...
            'linewidth',options.sim.lw,...
            'linestyle',options.sim.ls,...
            'color',options.sim.col); hold on;
        xlabel('time'); ylabel('mean');
        xlim([Data.SCTLstat.time(1),Data.SCTLstat.time(end)]);
    end

    % Covariances - Data
    I = [];
    for i = 1:ny
        I = [I,(2*i*ny+1):((2*i+1)*ny)];
    end
    Cz_abs_max = max([max(abs(Data.SCTLstat.Cz(:))),max(abs(Sim.Cz(:)))]);

    subplot(ny+1,2*ny,I); hold off;
    imagesc(Data.SCTLstat.Cz,Cz_abs_max*[-1,1]); hold on;
    xlabel('column'); ylabel('row');
    colorbar;
    caxis([min(min(min(Data.SCTLstat.Cz,Cz_abs_max))),max(max(max(Data.SCTLstat.Cz,Cz_abs_max)))])
    % Covariances - Simulation
    I = [];
    for i = 1:ny
        I = [I,((2*i+1)*ny+1):((2*i+2)*ny)];
    end
    subplot(ny+1,2*ny,I); hold off;
    imagesc(Sim.Cz,Cz_abs_max*[-1,1]); hold on;
    xlabel('column'); ylabel('row');
    colorbar;
    caxis([min(min(min(Data.SCTLstat.Cz,Cz_abs_max))),max(max(max(Data.SCTLstat.Cz,Cz_abs_max)))])

end

%% Visualization: Data
if isempty(Sim)
    % Number of outputs
    ny = size(Data.SCTLstat.mz,2);
    
    % Means
    for i = 1:ny
        subplot(ny+1,ny,i); hold off;
        plot(Data.SCTLstat.time,Data.SCTLstat.mz(:,i),'-',...
            'linewidth',options.data.lw,...
            'linestyle',options.data.ls,...
            'color',options.data.col); hold on;
        plot(Data.SCTLstat.time,Data.SCTLstat.mz(:,i)-3*Data.SCTLstat.Sigma_mz(:,i),'-',...
            'linewidth',options.error.lw,...
            'linestyle',options.error.ls,...
            'color',options.error.col); hold on;
        plot(Data.SCTLstat.time,Data.SCTLstat.mz(:,i)+3*Data.SCTLstat.Sigma_mz(:,i),'-',...
            'linewidth',options.error.lw,...
            'linestyle',options.error.ls,...
            'color',options.error.col); hold on;
        xlabel('time'); ylabel('mean');
        xlim([Data.SCTLstat.time(1),Data.SCTLstat.time(end)]);
    end

    % Covariances
    Cz_abs_max = max(abs(Data.SCTLstat.Cz(:)));
    subplot(ny+1,ny,(ny+1):(ny+1)*ny); hold off;
    imagesc(Data.SCTLstat.Cz,Cz_abs_max*[-1,1]); hold on;
    xlabel('column'); ylabel('row');
    colorbar;
end

%%
drawnow

end

