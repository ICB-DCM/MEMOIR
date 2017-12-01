% function plotPA(Data,Sim,s,options)
function plotPA(varargin)

    persistent fh

    %% Check and assign inputs
    if nargin >= 2
        Data = varargin{1};
        Sim = varargin{2};  

        % getSimulationMEMOIR passes also mFine, logL_PA only passes .m
        if ~isfield(Sim, 'mFine')
            Sim.mFine = Sim.m;
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
    options.data.ls = 'none';
    options.data.mean_lw = 2;
    options.data.bound_lw = 1;
    options.data.bound_ls = '--';
    options.sim.col = 'r';
    options.sim.area_col = [0.7,0.7,0.7];
    options.sim.border_col = [0.6,0.6,0.6];
    options.sim.ls = '-';
    options.sim.mean_lw = 2;
    options.sim.bound_lw = 1;
    options.sim_true.col = 'g';
    options.sim_true.area_col = [0.2,0.7,0.2];
    options.sim_true.border_col = [0,0.6,0];
    options.sim_true.ls = '--';
    options.sim_true.mean_lw = 2;
    options.sim_true.bound_lw = 1;
    options.error.col = 'b';
    options.error.ls = 'none';
    options.error.lw = 1;
    options.title = '';
    if nargin == 4
        options = setdefault(varargin{4},options);
    end
    set(gcf,'Name',options.title);

    %% Subplot dimensions
    n_y = size(Sim.m, 2);
    if ~isempty(Sim)
        nc = 2;
        nr = n_y;
    else
        nc = ceil(sqrt(n_y));
        nr = ceil(n_y/nc);
    end

    % Check if this experiment is a dose response experiment
    if (isfield(Data, 'condition') && ~isempty(Data.condition))
        conditions = unique(Data.condition,'rows');
        numConditions = size(conditions,1);
        if (numConditions > 1)
            % Get the label for the condition
            if (isfield(Data, 'conditionName') && ~isempty(Data.conditionName))
                conditionName = Data.conditionName;
            else
                warning('It looks like a dose-response-experiment should be plotted, but Data.conditionName is not set.');
                conditionName = 'dose';
            end

            % Find the row, which contains the information about the dose
            drRow = find(conditions(2,:) - conditions(1,:),1);
            doses = conditions(:,drRow);
            steps = Data.condition(:,drRow);
            tmpSteps = nan(size(steps));
            for iStep = 1:length(steps)
                for iDose = 1:length(doses)
                    if (steps(iStep) == doses(iDose))
                        tmpSteps(iStep) = iDose;
                    end
                end
            end
            dataSteps = tmpSteps;
            steps = (1 : length(doses))';
        else
            steps = Sim.t;
            dataSteps = Data.PA.time;
        end
    else
        numConditions = 1;
        if ~isempty(Sim)
            steps = Sim.t;
        end
        dataSteps = Data.PA.time;
    end

    %% Visualization: Data and Simulation
    if ~isempty(Sim)
        % Loop: measurands
        for j = 1:n_y
            % Data and simulation 
            subplot(nr, nc, 2*(j-1)+1); 
            hold off;

            lhCount = 1;

            % Plot noise
            lh(lhCount) = fill([steps(1:end); steps(end:-1:1)],...
                 [Sim.mFine(1:end,j) - Sim.Sigma_m(1:end,j);...
                  Sim.mFine(end:-1:1,j) + Sim.Sigma_m(end:-1:1,j)],...
                  options.sim.area_col);  
            hold on;
            plot(steps, Sim.mFine(:,j) - Sim.Sigma_m(:,j),'-',...
                  'linewidth',options.sim.bound_lw,...
                  'linestyle',options.sim.ls,...
                  'color',options.sim.border_col);
            plot(steps, Sim.mFine(:,j) +Sim.Sigma_m(:,j),'-',...
                  'linewidth',options.sim.bound_lw,...
                  'linestyle',options.sim.ls,...
                  'color',options.sim.border_col);
            alpha(0.5);

            if samples_vs_sp
                % Plot noise
                lhCount = lhCount + 1;
                lh(lhCount) = fill([steps(1:end); steps(end:-1:1)],...
                     [Sim.mFineTrue(1:end,j) - Sim.Sigma_m(1:end,j);...
                      Sim.mFineTrue(end:-1:1,j) + Sim.Sigma_m(end:-1:1,j)],...
                      options.sim_true.area_col);  
                plot(steps, Sim.mFineTrue(:,j) - Sim.Sigma_m(:,j),'-',...
                      'linewidth',options.sim_true.bound_lw,...
                      'linestyle',options.sim_true.ls,...
                      'color',options.sim_true.border_col);
                plot(steps, Sim.mFineTrue(:,j) + Sim.Sigma_m(:,j),'-',...
                      'linewidth',options.sim_true.bound_lw,...
                      'linestyle',options.sim_true.ls,...
                      'color',options.sim_true.border_col);
                alpha(0.5);
            end


            % Plot data points
            lhCount = lhCount + 1;
            lh(lhCount) = plot(dataSteps, Data.PA.m(:,j), '+',...
                  'linewidth',options.data.mean_lw,...
                  'linestyle',options.data.ls,...
                  'color',options.data.col);

            % Plot simulation
            lhCount = lhCount + 1;
            lh(lhCount) = plot(steps, Sim.mFine(:,j),'-',...
                  'linewidth',options.sim.mean_lw,...
                  'linestyle',options.sim.ls,...
                  'color',options.sim.col);

            if samples_vs_sp
                % Plot noise
                lhCount = lhCount + 1;
                lh(lhCount) = plot(steps, Sim.mFineTrue(:,j),'-',...
                    'linewidth',options.sim_true.mean_lw,...
                    'linestyle',options.sim_true.ls,...
                    'color',options.sim_true.col);
            end

            if (numConditions == 1)
                xlabel('time'); 
            else
                xlabel(conditionName);
                thisPlot = gca();
                thisPlot.XTick = 1:length(doses);
                thisPlot.XTickLabel = num2str(doses);
            end
            ylabel(Data.measurands{j});
            xlim(dataSteps([1,end]));
            if (j == 1)
                if samples_vs_sp
                    legend(lh,{'noise (SP)', 'noise (Sam)', 'data', 'model (SP)', 'model (Sam)'});
                else
                    legend(lh,{'noise', 'data', 'model'});
                end
            end

            subplot(nr,nc,2*(j-1)+2); 
            hold off;
            plot(dataSteps,Data.PA.m(:,j)-Sim.m(:,j),'o',...
                'linewidth',options.error.lw,...
                'linestyle',options.error.ls,...
                'color',options.error.col); hold on;

            if (numConditions == 1)
                xlabel('time'); 
            else
                xlabel(conditionName);
                thisPlot = gca();
                thisPlot.XTick = 1:length(doses);
                thisPlot.XTickLabel = num2str(doses);
            end
            ylabel(['error of mean(' Data.measurands{j} ')']);
            xlim(steps([1,end]));
        end
    end

    %% Visualization: Data
    if isempty(Sim)
        % Loop: measurands
        for j = 1:n_y
            subplot(nr,nc,j); 
            hold off;
            plot(dataSteps, Data.PA.m(:,j), 'o',...
                'linewidth',options.data.mean_lw,...
                'linestyle',options.data.ls,...
                'color',options.data.col); 
            hold on;

            if (numConditions == 1)
                xlabel('time'); 
            else
                xlabel(conditionName);
                thisPlot = gca();
                thisPlot.XTick = 1:length(doses);
                thisPlot.XTickLabel = num2str(doses);
            end
            ylabel(['error of var(' Data.measurands{j} ')']);
            xlim(dataSteps([1,end]));
        end
    end


    %%
    drawnow

end

