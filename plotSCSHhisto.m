% function fh = plotSCSH(Data,Sim,s,options)
function plotSCSHhisto(varargin)

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
    
    % If also sampling results should be plottet vs sigma points
    if ~isfield(Sim, 'mFineTrue')
        samples_vs_sp = 0;
    else
        samples_vs_sp = 1;
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
options.sim_true.col = 'g';
options.sim_true.area_col = [0.2,0.7,0.2];
options.sim_true.border_col = [0,0.6,0];
options.sim_true.ls = '--';
options.sim_true.mean_lw = 2;
options.sim_true.bound_lw = 1;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
options.title = '';
if nargin == 4
    options = setdefault(varargin{4},options);
end
set(gcf,'Name',options.title);

%% Subplot dimensions
n_ym = size(Data.SCSH.m, 2);
n_yC = size(Data.SCSH.C, 2);
time_vec = unique(Data.SCSH.time);
n_T = length(time_vec);
if ~isempty(Sim)
    nc = 3;
    nr = n_ym + n_yC;
else
    nc = ceil(sqrt(n_ym));
    nr = ceil(n_ym/nc);
end

%% Visualization: Data and Simulation
if ~isempty(Sim)
    lhCount = 0;
    % Loop: measurands - means
    for iTime = 1 : n_T
        for iMeas = 1 : n_ym
            % Data and simulation - mean
            subplot(n_T, n_ym, (iTime-1)*n_ym + iMeas); 
            hold off;

              
            lhCount = lhCount + 1;
            histogram(Sim.Y_true(iTime, iMeas, :), ...
                (0.1:0.0075:0.8), ...
                'Normalization', 'pdf');
            lh(lhCount) = gca;
            
            ylabel(['t = ' num2str(time_vec(iTime))]);
            if (iTime == n_T)
                xlabel(['Logicle of ' Data.measurands{iMeas}]);
            end
            if (iTime == 1)
                title([Data.measurands{iMeas} ' - Distribution']);
            end
            xlim([0.1, 0.8])
            ylim([0, 5])
        end
    end
end

nLH = lhCount;
for lhCount = 1 : nLH/2
    if lhCount < 10
       lh((lhCount-1)*2 + 1).XTick = [];
       lh((lhCount-1)*2 + 2).XTick = [];
    end
    lh((lhCount-1)*2 + 1).Position = [0.1, 0.95 - 0.09 * lhCount, 0.4, 0.085];
    lh((lhCount-1)*2 + 2).Position = [0.55, 0.95 - 0.09 * lhCount, 0.4, 0.085];
end
%         % Plot noise
%         lh(lhCount) = fill([Sim.t(1:end); Sim.t(end:-1:1)],...
%              [Sim.mFine(1:end,j) - Sim.Sigma_m(1:end,j);...
%               Sim.mFine(end:-1:1,j) + Sim.Sigma_m(end:-1:1,j)],...
%               options.sim_m.area_col);  
%         alpha(0.5);
%         hold on;
%         plot(Sim.t, Sim.mFine(:,j) - Sim.Sigma_m(:,j),'-',...
%               'linewidth',options.sim_m.bound_lw,...
%               'linestyle',options.sim_m.ls,...
%               'color',options.sim_m.border_col);
%         plot(Sim.t, Sim.mFine(:,j) + Sim.Sigma_m(:,j),'-',...
%               'linewidth',options.sim_m.bound_lw,...
%               'linestyle',options.sim_m.ls,...
%               'color',options.sim_m.border_col);
%           
%         if samples_vs_sp
%             % Plot noise
%             lhCount = lhCount + 1;
%             lh(lhCount) = fill([Sim.t(1:end); Sim.t(end:-1:1)],...
%                  [Sim.mFineTrue(1:end,j) - Sim.Sigma_m(1:end,j);...
%                   Sim.mFineTrue(end:-1:1,j) + Sim.Sigma_m(end:-1:1,j)],...
%                   options.sim_true.area_col);  
%             plot(Sim.t, Sim.mFineTrue(:,j) - Sim.Sigma_m(:,j),'-',...
%                   'linewidth',options.sim_true.bound_lw,...
%                   'linestyle',options.sim_true.ls,...
%                   'color',options.sim_true.border_col);
%             plot(Sim.t, Sim.mFineTrue(:,j) + Sim.Sigma_m(:,j),'-',...
%                   'linewidth',options.sim_true.bound_lw,...
%                   'linestyle',options.sim_true.ls,...
%                   'color',options.sim_true.border_col);
%             alpha(0.5);
%         end
%         
%         % Plot data
%         lhCount = lhCount + 1;
%         lh(lhCount) = plot(Data.SCSH.time, Data.SCSH.m(:,j),'+',...
%               'linewidth',options.data_m.mean_lw,...
%               'linestyle',options.data_m.ls,...
%               'color',options.data_m.col);
% 
%         % Plot simulation
%         lhCount = lhCount + 1;
%         lh(lhCount) = plot(Sim.t, Sim.mFine(:,j),'-',...
%               'linewidth',options.sim_m.mean_lw,...
%               'linestyle',options.sim_m.ls,...
%               'color',options.sim_m.col);
%         
%         if samples_vs_sp
%             % Plot simulation
%             lhCount = lhCount + 1;
%             lh(lhCount) = plot(Sim.t, Sim.mFineTrue(:,j),'-',...
%                 'linewidth',options.sim_true.mean_lw,...
%                 'linestyle',options.sim_true.ls,...
%                 'color',options.sim_true.col);
%         end
%           
%         xlabel('time'); ylabel(Data.measurands{j});
%         xlim(Data.SCSH.time([1,end]));
%         if (j == 1)
%             if samples_vs_sp
%                 legend(lh,{'noise (SP)', 'noise (Sam)', 'data', 'model mean (SP)', 'model meam (Sam)'});
%             else
%                 legend(lh,{'noise', 'data', 'model mean'});
%             end
%         end 
%         
%         % Error of mean
%         subplot(nr,nc,3*(j-1)+3); hold off;
%             plot(Data.SCSH.time,Data.SCSH.m(:,j)-Sim.m(:,j),'-',...
%                 'linewidth',options.error.lw,...
%                 'linestyle',options.error.ls,...
%                 'color',options.error.col); hold on;
%             
%         xlabel('time'); ylabel(['error of mean(' Data.measurands{j} ')']);
%         xlim(Data.SCSH.time([1,end]));
%     end
%     
%         % Data and simulation - mean
%         subplot(nr,nc,3*(n_ym + j-1)+[1,2]); 
%         hold off;
%         
%         % Plot noise of variability
%         lh2Count = 1;
%         lh2(lh2Count) = fill([Sim.t(1:end); Sim.t(end:-1:1)],...
%              [Sim.CFine(1:end,j) - Sim.Sigma_C(1:end,j);...
%               Sim.CFine(end:-1:1,j) + Sim.Sigma_C(end:-1:1,j)],...
%               options.sim_C.area_col); 
%         hold on;
%         alpha(0.5);
%         plot(Sim.t, Sim.CFine(:,j) - Sim.Sigma_C(:,j),'-',...
%               'linewidth',options.sim_C.bound_lw,...
%               'linestyle',options.sim_C.ls,...
%               'color',options.sim_C.border_col);
%         plot(Sim.t, Sim.CFine(:,j) + Sim.Sigma_C(:,j),'-',...
%               'linewidth',options.sim_C.bound_lw,...
%               'linestyle',options.sim_C.ls,...
%               'color',options.sim_C.border_col);
%           
%         if samples_vs_sp
%             % Plot noise
%             lh2Count = lh2Count + 1;
%             lh2(lh2Count) = fill([Sim.t(1:end); Sim.t(end:-1:1)],...
%                  [Sim.CFineTrue(1:end,j) - Sim.Sigma_C(1:end,j);...
%                   Sim.CFineTrue(end:-1:1,j) + Sim.Sigma_C(end:-1:1,j)],...
%                   options.sim_true.area_col);  
%             plot(Sim.t, Sim.CFineTrue(:,j) - Sim.Sigma_C(:,j),'-',...
%                   'linewidth',options.sim_true.bound_lw,...
%                   'linestyle',options.sim_true.ls,...
%                   'color',options.sim_true.border_col);
%             plot(Sim.t, Sim.CFineTrue(:,j) + Sim.Sigma_C(:,j),'-',...
%                   'linewidth',options.sim_true.bound_lw,...
%                   'linestyle',options.sim_true.ls,...
%                   'color',options.sim_true.border_col);
%             alpha(0.5);
%         end
%           
%         % Plot data 
%         lh2Count = lh2Count + 1;
%         lh2(lh2Count) = plot(Data.SCSH.time, Data.SCSH.C(:,j),'k+',...
%               'linewidth',options.data_C.mean_lw,...
%               'linestyle',options.data_C.ls,...
%               'color',options.data_C.col);
% 
%         % Plot simulation
%         lh2Count = lh2Count + 1;
%         lh2(lh2Count) = plot(Sim.t, Sim.CFine(:,j),'-',...
%               'linewidth',options.sim_C.mean_lw,...
%               'linestyle',options.sim_C.ls,...
%               'color',options.sim_C.col);
%           
%         if samples_vs_sp
%             % Plot simulation
%             lh2Count = lh2Count + 1;
%             lh2(lh2Count) = plot(Sim.t, Sim.CFineTrue(:,j),'-',...
%                 'linewidth',options.sim_true.mean_lw,...
%                 'linestyle',options.sim_true.ls,...
%                 'color',options.sim_true.col);
%         end
%           
%         xlabel('time'); ylabel([Data.measurands{j} ' - variability']);
%         xlim(Data.SCSH.time([1,end]));
%         if (j == 1)
%             if samples_vs_sp
%                 legend(lh2,{'noise (SP)', 'noise (Sam)', 'data variance', 'model variability (SP)', 'model variability (Sam)'});
%             else
%                 legend(lh2,{'noise', 'data variance', 'model variability'});
%             end
%         end 
%         
%         % Error of variance
%         subplot(nr,nc,3*(n_ym+j-1)+3); hold off;
%             plot(Data.SCSH.time,Data.SCSH.C(:,j)-Sim.C(:,j),'-',...
%                 'linewidth',options.error.lw,...
%                 'linestyle',options.error.ls,...
%                 'color',options.error.col); hold on;
%             
%         xlabel('time'); ylabel(['error of var(' Data.measurands{j} ')']);
%         xlim(Data.SCSH.time([1,end]));


%% Visualization: Data
if isempty(Sim)
    % Loop: measurands
    for j = 1:n_ym
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

