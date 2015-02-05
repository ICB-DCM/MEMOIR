% function fh = plotSCTL(Data,Sim,fh,options)
function fh = plotSCTL(varargin)

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
n_y = size(Data.SCTL.Y,2);
if ~isempty(Sim)
    nc = 2;
    nr = n_y;
else
    nc = ceil(sqrt(n_y));
    nr = ceil(n_y/nc);
end

%% Visualization: Data and Simulation of individuals
if ~isempty(Sim)
    % Loop: measurands
    for j = 1:n_y
        % Data and simulation 
        subplot(nr,nc,2*j-1); hold off;
        for i = 1:size(Data.SCTL.Y,3)
            lh(1) = plot(Data.SCTL.time,Data.SCTL.Y(:,j,i),'-',...
                'linewidth',options.data.lw,...
                'linestyle',options.data.ls,...
                'color',options.data.col); hold on;
            lh(2) = plot(Data.SCTL.time,Sim(:,j,i),'-',...
                'linewidth',options.sim.lw,...
                'linestyle',options.sim.ls,...
                'color',options.sim.col);
        end
        xlabel('time'); ylabel(Data.measurands{j});
        xlim(Data.SCTL.time([1,end]));
        if j == 1
            legend(lh,{'data','model'});
        end
        
        % Error 
        subplot(nr,nc,2*j); hold off;
        for i = 1:size(Data.SCTL.Y,3)
            plot(Data.SCTL.time,Data.SCTL.Y(:,j,i)-Sim(:,j,i),'-',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;
        end
        xlabel('time'); ylabel(['error ' Data.measurands{j}]);
        xlim(Data.SCTL.time([1,end]));
    end
end

%% Visualization: Data of individuals
if isempty(Sim)
    % Loop: measurands
    for j = 1:n_y
        % Data
        subplot(nr,nc,j); hold off;
        for i = 1:size(Data.SCTL.Y,3)
            plot(Data.SCTL.time,Data.SCTL.Y(:,j,i),'-',...
                'linewidth',options.data.lw,...
                'linestyle',options.data.ls,...
                'color',options.data.col); hold on;
        end
        xlabel('time'); ylabel(Data.measurands{j});
        xlim(Data.SCTL.time([1,end]));
    end
end

%%
drawnow

end

