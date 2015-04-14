% function fh = plotPA(Data,Sim,fh,options)
function fh = plotPAt0(varargin)

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
n_y = size(Data.PA.m,2);
if ~isempty(Sim)
    nc = 2;
    nr = n_y;
else
    nc = ceil(sqrt(n_y));
    nr = ceil(n_y/nc);
end

%% Visualization: Data and Simulation
if ~isempty(Sim)
    % Loop: measurands
    if max(Data.PA.time)==0 % only initial concentrations
        data = [Data.PA.m;Sim.m]';
        subplot(2,1,1)
        bar(data,'grouped');
        yl = [0 1.1*max(data(:))];
        set(gca,'Position',[.18 .6 .7 .3],'YLim',yl)
        set(gca,'XTickLabel',[]);
        xt=get(gca,'XTick');
        text(xt,-.05*yl(1,2)*ones(length(xt),1),Data.measurands,...
            'HorizontalAlignment','right','rotation',45);
        legend({'data','model'},'Location','NorthWest')
        subplot(2,1,2);
        data = Data.PA.m-Sim.m;
        bar(data);
        set(gca,'XTickLabel',[],'Position',[.18 .1 .7 .3],...
            'YLim',[-1.1 1.1]*max(abs(data(:))))
        legend({'error of PA quant.'},'Location','NorthWest')
    else
        for j = 1:n_y
            % Data and simulation 
            subplot(nr,nc,2*(j-1)+1); hold off;

            lh(1) = plot(Data.PA.time,Data.PA.m(:,j),'-',...
                  'linewidth',options.data.mean_lw,...
                  'linestyle',options.data.ls,...
                  'color',options.data.col); hold on;

            lh(2) = plot(Data.PA.time,Sim.m(:,j),'-',...
                  'linewidth',options.sim.mean_lw,...
                  'linestyle',options.sim.ls,...
                  'color',options.sim.col); hold on;

            xlabel('time'); ylabel(Data.measurands{j});
            xlim(Data.PA.time([1,end]));
            if j == 1
                legend(lh,{'data','model'});
            end

            subplot(nr,nc,2*(j-1)+2); hold off;
                plot(Data.PA.time,Data.PA.m(:,j)-Sim.m(:,j),'-',...
                    'linewidth',options.error.lw,...
                    'linestyle',options.error.ls,...
                    'color',options.error.col); hold on;

            xlabel('time'); ylabel(['error of mean(' Data.measurands{j} ')']);
            xlim(Data.PA.time([1,end]));
        end
    end
end

%% Visualization: Data
if isempty(Sim)
    % Loop: measurands
    if max(Data.PA.time)==0 % only initial concentrations
        bar(Data.PA.m);
        set(gca,'Position',[.18 .2 .7 .7])
        set(gca,'XTickLabel',[]);
        xt=get(gca,'XTick');
        yl = get(gca,'YLim');
        text(xt,ones(length(xt),1)*(yl(1)-.05*diff(yl)),Data.measurands,...
            'HorizontalAlignment','right','rotation',45);
    else
        for j = 1:n_y
            subplot(nr,nc,j); hold off;
            plot(Data.PA.time,Data.PA.m(:,j),'-',...
                  'linewidth',options.data.mean_lw,...
                  'linestyle',options.data.ls,...
                  'color',options.data.col); hold on;

            xlabel('time'); ylabel(['error of var(' Data.measurands{j} ')']);
            xlim(Data.PA.time([1,end]));
        end
    end
end

%%
drawnow

end

