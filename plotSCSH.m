% function fh = plotSCSH(Data,Sim,s,options)
function plotSCSH(varargin)

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
options.data.area_col = [0.7,0.7,1];
options.data.ls = '-';
options.data.mean_lw = 2;
options.data.bound_lw = 1;
options.sim.col = 'r';
options.sim.ls = '--';
options.sim.mean_lw = 2;
options.sim.bound_lw = 1;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% Subplot dimensions
n_y = size(Data.SCSH.m,2);
if ~isempty(Sim)
    nc = 3;
    nr = 2*n_y;
else
    nc = ceil(sqrt(n_y));
    nr = ceil(n_y/nc);
end

%% Visualization: Data and Simulation
if ~isempty(Sim)
    % Loop: measurands
    for j = 1:n_y
        % Data and simulation 
        subplot(nr,nc,6*(j-1)+[1,2,4,5]); hold off;
        
        fill([Data.SCSH.time(1:end);Data.SCSH.time(end:-1:1)],...
             [Data.SCSH.m(1:end,j)-sqrt(Data.SCSH.C(1:end,j,j));...
              Data.SCSH.m(end:-1:1,j)+sqrt(Data.SCSH.C(end:-1:1,j,j))],options.data.area_col); hold on;
        lh(1) = plot(Data.SCSH.time,Data.SCSH.m(:,j),'-',...
              'linewidth',options.data.mean_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col); hold on;
        plot(Data.SCSH.time,Data.SCSH.m(:,j)-sqrt(Data.SCSH.C(:,j,j)),'-',...
              'linewidth',options.data.bound_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col); hold on;
        plot(Data.SCSH.time,Data.SCSH.m(:,j)+sqrt(Data.SCSH.C(:,j,j)),'-',...
              'linewidth',options.data.bound_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col); hold on;

        lh(2) = plot(Data.SCSH.time,Sim.m(:,j),'-',...
              'linewidth',options.sim.mean_lw,...
              'linestyle',options.sim.ls,...
              'color',options.sim.col); hold on;
        plot(Data.SCSH.time,Sim.m(:,j)-sqrt(Sim.C(:,j,j)),'-',...
              'linewidth',options.sim.bound_lw,...
              'linestyle',options.sim.ls,...
              'color',options.sim.col); hold on;
        plot(Data.SCSH.time,Sim.m(:,j)+sqrt(Sim.C(:,j,j)),'-',...
              'linewidth',options.sim.bound_lw,...
              'linestyle',options.sim.ls,...
              'color',options.sim.col); hold on;
        
        xlabel('time'); ylabel(Data.measurands{j});
        xlim(Data.SCSH.time([1,end]));
        if j == 1
            legend(lh,{'data','model'});
        end
        
        % Error of mean
        subplot(nr,nc,6*(j-1)+3); hold off;
            plot(Data.SCSH.time,Data.SCSH.m(:,j)-Sim.m(:,j),'-',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;
            
        xlabel('time'); ylabel(['error of mean(' Data.measurands{j} ')']);
        xlim(Data.SCSH.time([1,end]));

        % Error of variance
        subplot(nr,nc,6*(j-1)+6); hold off;
            plot(Data.SCSH.time,Data.SCSH.C(:,j,j)-Sim.C(:,j,j),'-',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;
            
        xlabel('time'); ylabel(['error of var(' Data.measurands{j} ')']);
        xlim(Data.SCSH.time([1,end]));
    end
end

%% Visualization: Data
if isempty(Sim)
    % Loop: measurands
    for j = 1:n_y
        subplot(nr,nc,j); hold off;
        fill([Data.SCSH.time(1:end);Data.SCSH.time(end:-1:1)],...
             [Data.SCSH.m(1:end,j)-sqrt(Data.SCSH.C(1:end,j,j));...
              Data.SCSH.m(end:-1:1,j)+sqrt(Data.SCSH.C(end:-1:1,j,j))],options.data.area_col); hold on;
        lh(1) = plot(Data.SCSH.time,Data.SCSH.m(:,j),'-',...
              'linewidth',options.data.mean_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col); hold on;
        plot(Data.SCSH.time,Data.SCSH.m(:,j)-sqrt(Data.SCSH.C(1:end,j,j)),'-',...
              'linewidth',options.data.bound_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col); hold on;
        plot(Data.SCSH.time,Data.SCSH.m(:,j)+sqrt(Data.SCSH.C(1:end,j,j)),'-',...
              'linewidth',options.data.bound_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col); hold on;

        xlabel('time'); ylabel(['error of var(' Data.measurands{j} ')']);
        xlim(Data.SCSH.time([1,end]));
    end
end

%%
drawnow

end

