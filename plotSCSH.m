% function fh = plotSCSH(Data,Sim,s,options)
function plotSCSH(varargin)

persistent fh

%% Check and assign inputs
if nargin >= 2
    Data = varargin{1};
    Sim = varargin{2};
    
    % getSimulationMEMOIR passes also mFine, logL_PA only passes .m
    if ~isfield(Sim, 'mFine')
        Sim.mFine = Sim.m;
    end
    if ~isfield(Sim, 'CFine')
        Sim.CFine = Sim.C;
    end
else
    error('Not enough inputs.')
end

% Figure handle
s = varargin{3};
fh = [];
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
options.data_m.col = 'k';
options.data_m.ls = 'none';
options.data_m.mean_lw = 2;
options.data_C.col = 'k';
options.data_C.ls = 'none';
options.data_C.mean_lw = 2;
options.sim_m.col = 'r';
options.sim_m.ls = '-';
options.sim_m.border_col = [0.6,0.6,0.6];
options.sim_m.area_col = [0.7,0.7,0.7];
options.sim_m.mean_lw = 2;
options.sim_m.bound_lw = 1;
options.sim_C.col = [0,0,0.7];
options.sim_C.ls = '-';
options.sim_C.border_col = [0.3,0.4,0.9];
options.sim_C.area_col = [0.5,0.6,1];
options.sim_C.mean_lw = 2;
options.sim_C.bound_lw = 1;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
options.title = '';
if nargin == 4
    options = setdefault(varargin{4},options);
end
set(gcf,'Name',options.title);

%% Subplot dimensions
n_y = size(Data.SCSH.m, 2);
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
        % Data and simulation - mean
        subplot(nr,nc,6*(j-1)+[1,2]); hold off;
        
        % Plot noise
        lh(1) = fill([Sim.t(1:end); Sim.t(end:-1:1)],...
             [Sim.mFine(1:end,j) - Sim.Sigma_m(1:end,j);...
              Sim.mFine(end:-1:1,j) + Sim.Sigma_m(end:-1:1,j)],...
              options.sim_m.area_col);  
        alpha(0.5);
        hold on;
        plot(Sim.t, Sim.mFine(:,j) - Sim.Sigma_m(:,j),'-',...
              'linewidth',options.sim_m.bound_lw,...
              'linestyle',options.sim_m.ls,...
              'color',options.sim_m.border_col);
        plot(Sim.t, Sim.mFine(:,j) + Sim.Sigma_m(:,j),'-',...
              'linewidth',options.sim_m.bound_lw,...
              'linestyle',options.sim_m.ls,...
              'color',options.sim_m.border_col);
        
        % Plot data 
        lh(2) = plot(Data.SCSH.time, Data.SCSH.m(:,j),'+',...
              'linewidth',options.data_m.mean_lw,...
              'linestyle',options.data_m.ls,...
              'color',options.data_m.col);

        % Plot simulation
        lh(3) = plot(Sim.t, Sim.mFine(:,j),'-',...
              'linewidth',options.sim_m.mean_lw,...
              'linestyle',options.sim_m.ls,...
              'color',options.sim_m.col);
        
        xlabel('time'); ylabel(Data.measurands{j});
        xlim(Data.SCSH.time([1,end]));
        if j == 1
            legend(lh,{'noise', 'data', 'model mean'});
        end
        
        % Error of mean
        subplot(nr,nc,6*(j-1)+3); hold off;
            plot(Data.SCSH.time,Data.SCSH.m(:,j)-Sim.m(:,j),'-',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;
            
        xlabel('time'); ylabel(['error of mean(' Data.measurands{j} ')']);
        xlim(Data.SCSH.time([1,end]));

        % Data and simulation - covariance
        subplot(nr,nc,6*(j-1)+[4,5]); hold off;
        
        % Plot noise of variability
        lh2(1) = fill([Sim.t(1:end); Sim.t(end:-1:1)],...
             [Sim.CFine(1:end,j) - Sim.Sigma_C(1:end,j);...
              Sim.CFine(end:-1:1,j) + Sim.Sigma_C(end:-1:1,j)],...
              options.sim_C.area_col); 
        hold on;
        alpha(0.5);
        plot(Sim.t, Sim.CFine(:,j) - Sim.Sigma_C(:,j),'-',...
              'linewidth',options.sim_C.bound_lw,...
              'linestyle',options.sim_C.ls,...
              'color',options.sim_C.border_col);
        plot(Sim.t, Sim.CFine(:,j) + Sim.Sigma_C(:,j),'-',...
              'linewidth',options.sim_C.bound_lw,...
              'linestyle',options.sim_C.ls,...
              'color',options.sim_C.border_col);
          
        % Plot data 
        lh2(2) = plot(Data.SCSH.time, Data.SCSH.C(:,j),'k+',...
              'linewidth',options.data_C.mean_lw,...
              'linestyle',options.data_C.ls,...
              'color',options.data_C.col);

        % Plot simulation
        lh2(3) = plot(Sim.t, Sim.CFine(:,j),'-',...
              'linewidth',options.sim_C.mean_lw,...
              'linestyle',options.sim_C.ls,...
              'color',options.sim_C.col);
          
        xlabel('time'); ylabel([Data.measurands{j} ' - variability']);
        xlim(Data.SCSH.time([1,end]));
        if j == 1
            legend(lh2,{'noise', 'variance', 'model variability'});
        end
        
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

