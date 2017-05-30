% function plotPA(Data,Sim,s,options)
function plotPA(varargin)

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

% Simulations for sigma points
if nargin >= 4
    SP = varargin{4};
else
    SP = [];
end

% Options
options.data.col = 'b';
options.data.area_col = [0.7,0.7,0.7];
options.data.ls = 'none';
options.data.mean_lw = 2;
options.data.bound_lw = 1;
options.data.bound_ls = '--';
options.sim.col = 'r';
options.sim.ls = '-';
options.sim.mean_lw = 1;
options.sim.bound_lw = 1;
options.error.col = 'b';
options.error.ls = 'none';
options.error.lw = 1;
options.title = '';
if nargin == 4
    options = setdefault(varargin{4},options);
end
set(gcf,'Name',options.title);

%% Subplot dimensions
mIndices = Data.PA.measuredIndices;
n_y = size(mIndices, 2);
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
    for j = 1:n_y
        % Data and simulation 
        subplot(nr, nc, 2*(j-1)+1); 
        hold off;
        
        lh(1) = jbfill(Data.PA.time', ...
            Sim.m(:,j)' + Data.PA.Sigma_m(:,j)', ...
            Sim.m(:,j)' - Data.PA.Sigma_m(:,j)', ...
            options.data.area_col, ...
            options.data.area_col, ...
            0, 0.8);
        hold on;
         
        lh(2) = plot(Data.PA.time, Data.PA.m(:,mIndices(j)), '+',...
              'linewidth',options.data.mean_lw,...
              'linestyle',options.data.ls,...
              'color',options.data.col);

        lh(3) = plot(Data.PA.time, Sim.m(:,j),'-',...
              'linewidth',options.sim.mean_lw,...
              'linestyle',options.sim.ls,...
              'color',options.sim.col);
        
        xlabel('time'); ylabel(Data.measurands{j});
        xlim(Data.PA.time([1,end]));
        if j == 1
            legend(lh,{'data','model'});
        end
        
        subplot(nr,nc,2*(j-1)+2); hold off;
            plot(Data.PA.time,Data.PA.m(:,mIndices(j))-Sim.m(:,j),'o',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;
            
        xlabel('time'); ylabel(['error of mean(' Data.measurands{j} ')']);
        xlim(Data.PA.time([1,end]));
    end
end

%% Visualization: Data
if isempty(Sim)
    % Loop: measurands
    for j = 1:n_y
        subplot(nr,nc,j); 
        hold off;
        plot(Data.PA.time, Data.PA.m(:,mIndices(j)), 'o',...
            'linewidth',options.data.mean_lw,...
            'linestyle',options.data.ls,...
            'color',options.data.col); 
        hold on;

        xlabel('time'); ylabel(['error of var(' Data.measurands{j} ')']);
        xlim(Data.PA.time([1,end]));
    end
end

%%
drawnow

end

