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
if(isfield(Data.SCTL,'Y'))
    n_y = size(Data.SCTL.Y,2);
else
    n_y = 0;
end

if(isfield(Data.SCTL,'T'))
    n_t = size(Data.SCTL.T,1);
else
    n_t = 0;
end

n_o = n_y + n_t;

if ~isempty(Sim)
    nc = 2;
    nr = n_o;
else
    nc = ceil(sqrt(n_o));
    nr = ceil(n_o/nc);
end

%% Visualization: Data and Simulation of individuals
if ~isempty(Sim)
    % Loop: measurands
    for j = 1:n_y
        % Data and simulation 
        subplot(nr,nc,2*j-1); hold off;
        for i = 1:size(Data.SCTL.Y,3)
            ind = ~isnan(Data.SCTL.Y(:,j,i));
            lh(1) = plot(Data.SCTL.time(ind,1),Data.SCTL.Y(ind,j,i),'-',...
                'linewidth',options.data.lw,...
                'linestyle',options.data.ls,...
                'color',options.data.col); hold on;
            lh(2) = plot(Data.SCTL.time,Sim.Y(:,j,i),'-',...
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
            ind = ~isnan(Data.SCTL.Y(:,j,i));
            plot(Data.SCTL.time(ind,1),Data.SCTL.Y(ind,j,i)-Sim(ind,j,i),'-',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;
        end
        xlabel('time'); ylabel(['error ' Data.measurands{j}]);
        xlim(Data.SCTL.time([1,end]));
    end
    % Loop: events
    for j = 1:n_t
        % Data and simulation 
        subplot(nr,nc,2*(n_y+j)-1); hold off;
        for i = 1:size(Data.SCTL.T,3)
            ind = ~isnan(Data.SCTL.T(:,j,i));
            lh(1) = stem(Data.SCTL.T(ind,j,i),'-',...
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
            ind = ~isnan(Data.SCTL.Y(:,j,i));
            plot(Data.SCTL.time(ind,1),Data.SCTL.Y(ind,j,i)-Sim(ind,j,i),'-',...
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
            ind = ~isnan(Data.SCTL.Y(:,j,i));
            plot(Data.SCTL.time(ind,1),Data.SCTL.Y(ind,j,i),'-',...
                'linewidth',options.data.lw,...
                'linestyle',options.data.ls,...
                'color',options.data.col); hold on;
        end
        xlabel('time'); ylabel(Data.measurands{j});
        xlim(Data.SCTL.time([1,end]));
    end
    % Loop: events
    for j = 1:n_t
        % Data
        subplot(nr,nc,n_y+j); hold off;
        for i = 1:size(Data.SCTL.T,3)
            ind = ~isnan(Data.SCTL.T(j,:,i));
            stem(Data.SCTL.T(j,ind,i),zeros(size(Data.SCTL.T(j,ind,i))),...
                'color',options.data.col); hold on;
        end
        xlabel('time'); ylabel(Data.events{j});
        xlim(Data.SCTL.time([1,end]));
    end
end

%%
drawnow

end

