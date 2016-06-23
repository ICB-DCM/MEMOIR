function [P,logL_sc,dlogL_scdxi,ddlogL_scdxidxi] = logL_SCTL(xi, model, data, s, options, P)

persistent fp
persistent fl

options.nderiv = nargout-2;

%% Construct fixed effects and covariance matrix
beta = model.beta(xi);
delta = model.delta(xi);

[D,~,~,~,~,~] = xi2D(delta,options.type_D);
% debugging:
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,1,3)
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,3,5)
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,2,4)
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(delta,@(x) xi2D(x,options.type_D),1e-4,4,6)

% Initialization measurements
Sim_SCTL.Y = NaN(size(data.SCTL.Y));
% events
if(~isfield(data.SCTL,'T'))
    data.SCTL.T = zeros(0,1,size(data.SCTL.Y,3));
end
Sim_SCTL.T = NaN(size(data.SCTL.T));
Sim_SCTL.R = NaN(size(data.SCTL.T));

% set default scaling
if(~isfield(model,'SCTLscale'))
    model.SCTLscale = 1;
end

% Loop: Indiviudal cells

logLi_D = zeros(size(data.SCTL.Y,3),1);
logLi_T = zeros(size(data.SCTL.Y,3),1);
logLi_b = zeros(size(data.SCTL.Y,3),1);
logLi_I = zeros(size(data.SCTL.Y,3),1);

if options.nderiv >= 1
    dlogLi_Ddxi = zeros(size(data.SCTL.Y,3),length(xi));
    dlogLi_Tdxi = zeros(size(data.SCTL.Y,3),length(xi));
    dlogLi_bdxi = zeros(size(data.SCTL.Y,3),length(xi));
    dlogLi_Idxi = zeros(size(data.SCTL.Y,3),length(xi));
    if options.nderiv > 2
        ddlogLi_Ddxidxi = zeros(size(data.SCTL.Y,3),length(xi),length(xi));
        ddlogLi_Tdxidxi = zeros(size(data.SCTL.Y,3),length(xi),length(xi));
        ddlogLi_bdxidxi = zeros(size(data.SCTL.Y,3),length(xi),length(xi));
        ddlogLi_Idxidxi = zeros(size(data.SCTL.Y,3),length(xi),length(xi));
    end
end

tmp = arrayfun(@(x) any(~isnan(data.SCTL.Y(:,:,x)),2),1:size(data.SCTL.Y,3),'UniformOutput',false);
data.SCTL.ind_y = [tmp{:}];
tmp = arrayfun(@(x) any(~isnan(data.SCTL.T(:,:,x)),2),1:size(data.SCTL.T,3),'UniformOutput',false);
data.SCTL.ind_t = [tmp{:}];

Sim_SCTL_Y = zeros([length(data.SCTL.ind_y),size(data.SCTL.Y,2),size(data.SCTL.Y,3)]);
Sim_SCTL_T = zeros([length(data.SCTL.ind_t),size(data.SCTL.T,2),size(data.SCTL.T,3)]);
Sim_SCTL_R = zeros([length(data.SCTL.ind_t),size(data.SCTL.T,2),size(data.SCTL.T,3)]);
% check wether there is a parallel pool available
try
    p = gcp('nocreate')
catch
    p = [];
end
if isempty(p)
    for i = 1:size(data.SCTL.Y,3)
        [ logL_D,logL_T,logL_b,logL_I,bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i);
        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logL_D_SCTL_si(xi, model, data, s, options, P, i),1e-3,'val','dxi')
        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logL_T_SCTL_si(xi, model, data, s, options, P, i),1e-3,'val','dxi')
        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logL_b_SCTL_si(xi, model, data, s, options, P, i),1e-3,'val','dxi')
        % [g,g_fd_f,g_fd_b,g_fd_c]=testGradient(xi,@(xi) logL_I_SCTL_si(xi, model, data, s, options, P, i),1e-3,'val','dxi')
        logLi_D(i,1) = logL_D.val;
        logLi_T(i,1) = logL_T.val;
        logLi_b(i,1) = logL_b.val;
        if(options.integration)
            logLi_I(i,1) = logL_I.val;
        end
        b(:,i) = bhat.val;
        if(options.nderiv>0)
            dlogLi_Ddxi(i,:) = logL_D.dxi;
            dlogLi_Tdxi(i,:) = logL_T.dxi;
            dlogLi_bdxi(i,:) = logL_b.dxi;
            if(options.integration)
                dlogLi_Idxi(i,:) = logL_I.dxi;
            end
            dbdxi(i,:,:) = bhat.dxi;
            if(options.nderiv>1)
                ddlogLi_Ddxidxi(i,:,:) = logL_D.dxidxi;
                ddlogLi_Tdxidxi(i,:,:) = logL_T.dxidxi;
                ddlogLi_bdxidxi(i,:,:) = logL_b.dxidxi;
                if(options.integration)
                    ddlogLi_Idxidxi(i,:,:) = logL_I.dxidxi;
                end
            end
        end
        Sim_SCTL_Y(data.SCTL.ind_y(:,i),:,i) = Sim.SCTL_Y;
        Sim_SCTL_T(data.SCTL.ind_t(:,i),:,i) = Sim.SCTL_T;
        Sim_SCTL_R(data.SCTL.ind_t(:,i),:,i) = Sim.SCTL_R;
        
    end
else
    parfor i = 1:size(data.SCTL.Y,3)
        YY = zeros(size(data.SCTL.Y(:,:,i)))
        TT = zeros(size(data.SCTL.T(:,:,i)))
        RR = zeros(size(data.SCTL.T(:,:,i)))
        [ logL_D,logL_T,logL_b,logL_I,bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i)
        logLi_D(1,i) = logL_D.val;
        logLi_T(1,i) = logL_T.val;
        logLi_b(1,i) = logL_b.val;
        if(options.integration)
            logLi_I(1,i) = logL_I.val;
        end
        b(:,i) = bhat.val;
        if(options.nderiv>0)
            dlogLi_Ddxi(i,:) = logL_D.dxi;
            dlogLi_Tdxi(i,:) = logL_T.dxi;
            dlogLi_bdxi(i,:) = logL_b.dxi;
            if(options.integration)
                dlogLi_Idxi(i,:) = logL_I.dxi;
            end
            dbdxi(i,:,:) = bhat.dxi;
            if(options.nderiv>1)
                ddlogLi_Ddxidxi(i,:,:) = logL_D.dxidxi;
                ddlogLi_Tdxidxi(i,:,:) = logL_T.dxidxi;
                ddlogLi_bdxidxi(i,:,:) = logL_b.dxidxi;
                if(options.integration)
                    ddlogLi_Idxidxi(i,:,:) = logL_I.dxidxi;
                end
            end
        end
        YY(data.SCTL.ind_y(:,i),:) = Sim.SCTL_Y;
        TT(data.SCTL.ind_T(:,i),:) = Sim.SCTL_T;
        RR(data.SCTL.ind_T(:,i),:) = Sim.SCTL_R;
        Sim_SCTL_Y(:,:,i) = YY;
        Sim_SCTL_T(:,:,i) = TT;
        Sim_SCTL_R(:,:,i) = RR;
    end
end

P{s}.SCTL.bhat = b;
if(options.nderiv>0)
    P{s}.SCTL.dbhatdxi = dbdxi;
end
Sim_SCTL.Y = Sim_SCTL_Y;
Sim_SCTL.T = Sim_SCTL_T;
Sim_SCTL.R = Sim_SCTL_R;

%% Visulization
if options.plot
    
    % Visualisation of single cell parameters
    if(isempty(fp))
        if(isfield(model,'title'))
            if(ischar(model.title))
                fp(s) = figure('Name',model.title);
            else
                fp(s) = figure;
            end
        else
            fp(s) = figure;
        end
    else
        if(length(fp)<s)
            if(isfield(model,'title'))
                if(ischar(model.title))
                    fp(s) = figure('Name',model.title);
                else
                    fp(s) = figure;
                end
            else
                fp(s) = figure;
            end
        elseif(isempty(fp(s)))
            if(isfield(model,'title'))
                if(ischar(model.title))
                    fp(s) = figure('Name',model.title);
                else
                    fp(s) = figure;
                end
            else
                fp(s) = figure;
            end
        end
    end
    figure(fp(s))
    clf
    b_s = P{s}.SCTL.bhat;
    n_b = size(b_s,1);
    
    for j = 1:n_b
        subplot(ceil((n_b+1)/4),4,j+1)
        xx = linspace(-5*sqrt(D(j,j)),5*sqrt(D(j,j)),100);
        %nhist(P{s}.SCTL.bhat(j,:),'pdf','noerror');
        hold on
        plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','LineWidth',2)
        ecdf = zeros(length(xx),1);
        for k = 1:length(xx)
            ecdf(k) = sum(b_s(j,:)<xx(k))/length(b_s(j,:));
        end
        plot(xx,ecdf,'--r','LineWidth',2)
        
        
        if(j==1)
            
        end
        xlim([-5*sqrt(D(j,j)),5*sqrt(D(j,j))])
        ylim([0,1.1])
        %xlabel(char(model.sym.b(model.ind_b(j))));
        ylabel('cdf')
        box on
    end
    subplot(ceil(n_b+1/4),4,1,'Visible','off')
    hold on
    plot(xx,normcdf(xx,0,sqrt(D(j,j))),'.-b','Visible','off')
    plot(xx,ecdf,'--r','LineWidth',2,'Visible','off')
    
    
    legend('cdf of single cell Parameters','cdf of populaton Parameters')
    
    % Visualisation of likelihood contribution
    if(isempty(fl))
        if(isfield(model,'title'))
            if(ischar(model.title))
                fl(s) = figure('Name',model.title);
            else
                fl(s) = figure;
            end
        else
            fl(s) = figure;
        end
    else
        if(length(fl)<s)
            if(isfield(model,'title'))
                if(ischar(model.title))
                    fl(s) = figure('Name',model.title);
                else
                    fl(s) = figure;
                end
            else
                fl(s) = figure;
            end
        elseif(isempty(fl(s)))
            if(isfield(model,'title'))
                if(ischar(model.title))
                    fl(s) = figure('Name',model.title);
                else
                    fl(s) = figure;
                end
            else
                fl(s) = figure;
            end
        end
    end
    figure(fl(s))
    clf
    if(options.integration)
        bar(transpose([logLi_D,logLi_T,logLi_b,logLi_I]),'stacked')
        set(gca,'XTickLabel',{'data','event','par','int'})
    else
        bar(transpose([logLi_D,logLi_T,logLi_b]),'stacked')
        set(gca,'XTickLabel',{'data','event','par'})
    end
    ylabel('log-likelihood')
    title('likelihood contribution')
    
    
    
    % Visualisation of data and fit
    model.plot(data,Sim_SCTL,s);
end

%% Summation
logL_sc = logLi_D + logLi_b;
if(options.events)
    logL_sc = logL_sc + logLi_T;
end
if(options.integration)
    logL_sc = logL_sc + logLi_I;
end
if(options.nderiv >= 1)
    dlogL_scdxi = dlogLi_Ddxi + dlogLi_bdxi;
    if(options.events)
        dlogL_scdxi = dlogL_scdxi + dlogLi_Tdxi;
    end
    if(options.integration)
        dlogL_scdxi = dlogL_scdxi + dlogLi_Idxi;
    end
    
    if(options.nderiv >= 2)
        ddlogL_scdxidxi = ddlogLi_Ddxidxi + ddlogLi_bdxidxi;
        if(options.events)
            ddlogL_scdxidxi = ddlogL_scdxidxi + ddlogLi_Tdxidxi;
        end
        if(options.integration)
            ddlogL_scdxidxi = ddlogL_scdxidxi + ddlogLi_Idxidxi;
        end
    end
end


end

%%% FOR DEBUGGING

function logL_D = logL_D_SCTL_si(xi, model, data, s, options, P, i)

[ logL_D,logL_T,logL_b,logL_I,bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i);

end

function logL_T = logL_T_SCTL_si(xi, model, data, s, options, P, i)

[ logL_D,logL_T,logL_b,logL_I,bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i);

end

function logL_b = logL_b_SCTL_si(xi, model, data, s, options, P, i)

[ logL_D,logL_T,logL_b,logL_I,bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i);

end

function logL_I = logL_I_SCTL_si(xi, model, data, s, options, P, i)

[ logL_D,logL_T,logL_b,logL_I,bhat, Sim ] = logL_SCTL_si(xi, model, data, s, options, P, i);

end