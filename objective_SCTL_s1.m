% objective_SCTL_s1 is an auxiliary function for logL_CE_w_grad_2 and
% serves as objective function for the estimation of random effect
% parameters
%
% USAGE:
% ======
% [J,J.db,J.dbdb,...] = objective_SCTL_s1(model,beta,b,data.condition,delta,options.type_D,t,Ym,Tm,ind_y,ind_t,s)
%
% INPUTS:
% =======
%
% model ... model definition
% beta ... common effect parameter
% b ... random effect parameter
% data.condition ... experimental condition
% delta ... parametrisation of random effect covariance
% options.type_D ... covariance parametrisation type
% t ... time vector for simulation
% Ym ... measurements
% Tm ... observed event-times
% ind_y ... indexing of measurements
% ind_t ... indexing of events
% s ... experimental index
%
% Outputs:
% ========
% objective function J and derivatives wrt to b, beta and delta, see file
% for details, pd indicates that only the partial derivative is considered
% J
% J.db
% J.dbdb = G
% J.dbdbeta
% J.dbddelta
% J.dbetadbeta
% J.ddeltaddelta
% J.dbetaddelta
% dGdb
% pdGpdbeta
% pdGpddelta
% J.dbdbetadbeta
% J.dbddeltaddelta
% J.dbdbetaddelta
% ddGdbdb
% pddGdbpdbeta
% pdpdGpdbetapdbeta
% pddGdbpddelta
% pdpdGpddeltapddelta
% pdpdGpdbetapddelta
%
% 2015/04/14 Fabian Froehlich

function varargout = objective_SCTL_s1(model,data,beta,b,delta,s,i,options,nderiv)

[D,invD,dDddelta,dinvDddelta,ddDddeltaddelta,ddinvDddeltaddelta] = xi2D(delta,options.type_D);

% mixed effect parameter
phi.val = model.phi(beta,b);

nderiv = max(nderiv,nargout>1);

if(nargout>5)
    [J_D,J_T,Sim] = objective_phi(model,data,phi.val,s,i,options,nderiv,nargout>=3);
else
    [J_D,J_T] = objective_phi(model,data,phi.val,s,i,options,nderiv,nargout>=3);
end

% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(phi.val,@(phi) objective_phi_J_D(model,data,phi,s,i,options,nderiv),1e-6,'val','dphi')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(phi.val,@(phi) objective_phi_J_D(model,data,phi,s,i,options,nderiv),1e-6,'dphi','dphidphi')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(phi.val,@(phi) objective_phi_J_D(model,data,phi,s,i,options,nderiv),1e-6,'dphi','FIM')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(phi.val,@(phi) objective_phi_J_T(model,data,phi,s,i,options,nderiv),1e-6,'val','dphi')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(phi.val,@(phi) objective_phi_J_T(model,data,phi,s,i,options,nderiv),1e-6,'dphi','dphidphi')
% [g,g_fd_b,g_fd_f,g_fd_c] = testGradient(phi.val,@(phi) objective_phi_J_T(model,data,phi,s,i,options,nderiv),1e-6,'dphi','FIM')


switch(model.parameter_model)
    case 'normal'
        paramdist = @normal_param;
end
J_b = paramdist(b,delta,options.type_D,nderiv+(nargout>=3));

J.val = J_D.val + J_T.val + J_b.val ;
J.J_D = J_D;
J.J_T = J_T;
J.J_b = J_b;

if nargout >= 2
    
    phi.db = model.dphidb(beta,b);
    phi.dbeta = model.dphidbeta(beta,b);
    
    %% J.db
    
    J_D.db = chainrule(J_D.dphi,phi.db);
    
    J_T.db = chainrule(J_T.dphi,phi.db);
    
    J.db = J_D.db + J_T.db + J_b.db;
    
    %% J.dbeta
    
    J_D.dbeta = chainrule(J_D.dphi,phi.dbeta);
    
    J_T.dbeta = chainrule(J_T.dphi,phi.dbeta);
    
    J.dbeta = J_D.dbeta + J_T.dbeta;
    
    %% J.ddelta
    
    J.ddelta = J_b.ddelta;
    
    if nargout >= 3 || nderiv >= 1
        
        phi.dbdb = model.ddphidbdb(beta,b);
        %% J.dbdb
        % we need to make two different computations here,
        % one for the integration, in order to ensure that it is possible
        % to compute derivatives by using second order sensitivities by
        % using the FIM approximation and one that is used in the
        % computation of implicit derivatives where the accuracy is more
        % important
        
        J_D.appdbdphi = permute(sum(bsxfun(@times,J_D.FIM,permute(phi.db,[3,1,2])),2),[1,3,2]);
        J_D.appdbdb = permute(sum(bsxfun(@times,J_D.appdbdphi,permute(phi.db,[1,3,2])),1),[3,2,1]) ...
            + permute(chainrule(J_D.dphi,phi.dbdb),[2,3,1]);
        
        J_T.appdbdphi = permute(sum(bsxfun(@times,J_T.FIM,permute(phi.db,[3,1,2])),2),[1,3,2]);
        J_T.appdbdb = permute(sum(bsxfun(@times,J_T.appdbdphi,permute(phi.db,[1,3,2])),1),[3,2,1]) ...
            + permute(chainrule(J_T.dphi,phi.dbdb),[2,3,1]);
        
        % FIM
        FIM.val = J_D.appdbdb + J_T.appdbdb + J_b.dbdb;
        
        if nderiv >= 2
            
            phi.dbdbeta = model.ddphidbdbeta(beta,b);
            phi.dbetadbeta = model.ddphidbetadbeta(beta,b);
            
            J_D.dbdphi = permute(sum(bsxfun(@times,J_D.dphidphi,permute(phi.db,[3,1,2])),2),[3,1,2]);
            J_D.dbdb = permute(sum(bsxfun(@times,J_D.dbdphi,permute(phi.db,[3,1,2])),2),[1,3,2]) ...
                + permute(chainrule(J_D.dphi,phi.dbdb),[2,3,1]);
            
            J_T.dbdphi = permute(sum(bsxfun(@times,J_T.dphidphi,permute(phi.db,[3,1,2])),2),[3,1,2]);
            J_T.dbdb = permute(sum(bsxfun(@times,J_T.dbdphi,permute(phi.db,[3,1,2])),2),[1,3,2]) + permute(chainrule(J_T.dphi,phi.dbdb),[2,3,1]);
            
            J.dbdb = J_D.dbdb + J_T.dbdb + J_b.dbdb;
            
            %% J.dbdbeta
            
            J_D.dbdbeta = permute(sum(bsxfun(@times,J_D.dbdphi,permute(phi.dbeta,[3,1,2])),2),[1,3,2]) + permute(chainrule(J_D.dphi,phi.dbdbeta),[2,3,1]);
            J_T.dbdbeta = permute(sum(bsxfun(@times,J_T.dbdphi,permute(phi.dbeta,[3,1,2])),2),[1,3,2]) + permute(chainrule(J_T.dphi,phi.dbdbeta),[2,3,1]);
            
            J.dbdbeta = J_D.dbdbeta + J_T.dbdbeta;
            
            %% J.dbdelta
            
            J.dbddelta = J_b.dbddelta;
            
            %% J.dbetadbeta
            
            J_D.dbetadphi = permute(sum(bsxfun(@times,J_D.dphidphi,permute(phi.dbeta,[3,1,2])),2),[3,2,1]);
            J_T.dbetadphi = permute(sum(bsxfun(@times,J_T.dphidphi,permute(phi.dbeta,[3,1,2])),2),[3,2,1]);
            J_D.dbetadbeta = permute(sum(bsxfun(@times,J_D.dbetadphi,permute(phi.dbeta,[3,1,2])),2),[1,3,2]) + permute(chainrule(J_D.dphi,phi.dbetadbeta),[2,3,1]);
            J_T.dbetadbeta = permute(sum(bsxfun(@times,J_T.dbetadphi,permute(phi.dbeta,[3,1,2])),2),[1,3,2]) + permute(chainrule(J_T.dphi,phi.dbetadbeta),[2,3,1]);
            
            J.dbetadbeta = J_D.dbetadbeta + J_T.dbetadbeta;
            
            %% J.ddeltaddelta
            
            J.ddeltaddelta = J_b.ddeltaddelta;
            
            %% J.dbetaddelta
            
            J.dbetaddelta = zeros(length(beta),size(dDddelta,3));
            
            if nderiv >= 3
                
                %% FIM.db

                J_D.dbdphidphi = permute(sum(bsxfun(@times,J_D.FIMdphi,permute(phi.db,[1,3,4,2])),1),[4,2,3,1]);
                J_D.dbdbdphi = permute(sum(bsxfun(@times,J_D.dbdphidphi,permute(phi.db,[3,1,4,2])),2),[1,4,3,2]) ...
                    + permute(sum(bsxfun(@times,J_D.dphidphi,permute(phi.dbdb,[1,4,2,3])),1),[3,4,2,1]);
                J_D.dbdbdb = permute(sum(bsxfun(@times,J_D.dbdbdphi,permute(phi.db,[3,4,1,2])),3),[1,2,4,3]) ...
                    + permute(chainrule(J_D.dbdphi,phi.dbdb),[2,3,1]);
                
                J_T.dbdphidphi = permute(sum(bsxfun(@times,J_T.FIMdphi,permute(phi.db,[1,3,4,2])),1),[4,2,3,1]);
                J_T.dbdbdphi = permute(sum(bsxfun(@times,J_T.dbdphidphi,permute(phi.db,[3,1,4,2])),2),[1,4,3,2]) ...
                    + permute(sum(bsxfun(@times,J_T.dphidphi,permute(phi.dbdb,[1,4,2,3])),1),[3,4,2,1]);
                J_T.dbdbdb = permute(sum(bsxfun(@times,J_T.dbdbdphi,permute(phi.db,[3,4,1,2])),3),[1,2,4,3]) ...
                    + permute(chainrule(J_T.dbdphi,phi.dbdb),[2,3,1]);
                
                FIM.db = J_D.dbdbdb + J_T.dbdbdb + J_b.dbdbdb;
                
                %% FIM.dbeta
                
                J_D.dbdbdbeta = permute(sum(bsxfun(@times,J_D.dbdbdphi,permute(phi.dbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                J_T.dbdbdbeta = permute(sum(bsxfun(@times,J_T.dbdbdphi,permute(phi.dbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                FIM.dbeta = J_D.dbdbdbeta + J_T.dbdbdbeta;
                
                %% FIM.ddelta
                FIM.ddelta = J_b.dbdbddelta;
                
                %% J.dbdbetadbeta
                
                J_D.dbdbetadphi = permute(sum(bsxfun(@times,J_D.dbdphidphi,permute(phi.dbeta,[3,1,4,2])),2),[1,4,3,2]);
                J_D.dbdbetadbeta = permute(sum(bsxfun(@times,J_D.dbdbetadphi,permute(phi.dbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                J_T.dbdbetadphi = permute(sum(bsxfun(@times,J_T.dbdphidphi,permute(phi.dbeta,[3,1,4,2])),2),[1,4,3,2]);
                J_T.dbdbetadbeta = permute(sum(bsxfun(@times,J_T.dbdbetadphi,permute(phi.dbeta,[3,4,1,2])),3),[1,2,4,3]);
                
                J.dbdbetadbeta = J_D.dbdbetadbeta + J_T.dbdbetadbeta;
                
                %% J.dbddeltadelta
                J.dbddeltaddelta = J_b.dbddeltaddelta;
                
                %% J.dbdbetaddelta
                
                J.dbdbetaddelta = zeros(length(b),length(beta),size(dDddelta,3));
                
%                 if nderiv >= 4
%                     
%                     %                         %% FIM.dbdb
%                     %
%                     %                         temp = squeeze(sum(bsxfun(@times,bsxfun(@times,J_D.dYdY,permute(Y.dphidphi,[4,1,2,3])),permute(Y.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_D.dYdSigma,permute(Y.dphidphi,[4,1,2,3])),permute(Sigma_noise.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_D.dSigmadSigma,permute(Sigma_noise.dphidphi,[4,1,2,3])),permute(Sigma_noise.dphidphi,[4,1,5,6,2,3])),2));
%                     %
%                     %                         J_D..dphidphidphidphi = ...
%                     %                             permute(temp,[1,3,2,4]) + permute(temp,[1,3,4,2]) + permute(temp,[3,4,1,2]);
%                     %
%                     %                         J_D..dphidYdYdY = bsxfun(@times,J_D..dYdYdYdY,permute(Y.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_D..dYdYdYdSigma,permute(Sigma_noise.dphi,[3,1,2]));
%                     %                         J_D..dphidYdYdSigma = bsxfun(@times,J_D..dYdYdYdSigma,permute(Y.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_D..dYdYdSigmadSigma,permute(Sigma_noise.dphi,[3,1,2]));
%                     %                         J_D..dphidYdSigmadSigma = bsxfun(@times,J_D..dYdYdSigmadSigma,permute(Y.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_D..dYdSigmadSigmadSigma,permute(Sigma_noise.dphi,[3,1,2]));
%                     %                         J_D..dphidSigmadSigmadSigma = bsxfun(@times,J_D..dYdSigmadSigmadSigma,permute(Y.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_D..dSigmadSigmadSigmadSigma,permute(Sigma_noise.dphi,[3,1,2]));
%                     %                         J_D..dphidphidYdY = bsxfun(@times,J_D..dphidYdYdY,permute(Y.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_D..dphidYdYdSigma,permute(Sigma_noise.dphi,[3,1,4,2]));
%                     %                         J_D..dphidphidYdSigma = bsxfun(@times,J_D..dphidYdYdSigma,permute(Y.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_D..dphidYdSigmadSigma,permute(Sigma_noise.dphi,[3,1,4,2]));
%                     %                         J_D..dphidphidSigmadSigma = bsxfun(@times,J_D..dphidYdSigmadSigma,permute(Y.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_D..dphidSigmadSigmadSigma,permute(Sigma_noise.dphi,[3,1,4,2]));
%                     %                         J_D..dphidphidphidY = bsxfun(@times,J_D..dphidphidYdY,permute(Y.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_D..dphidphidYdSigma,permute(Sigma_noise.dphi,[3,1,4,5,2]));
%                     %                         J_D..dphidphidphidSigma = bsxfun(@times,J_D..dphidphidYdSigma,permute(Y.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_D..dphidphidSigmadSigma,permute(Sigma_noise.dphi,[3,1,4,5,2]));
%                     %                         J_D..dphidphidphidphi = J_D..dphidphidphidphi + squeeze(sum(bsxfun(@times,J_D..dphidphidphidY,permute(Y.dphi,[3,1,4,5,6,2])) ...
%                     %                             + bsxfun(@times,J_D..dphidphidphidSigma,permute(Sigma_noise.dphi,[3,1,4,5,6,2])),2));
%                     %
%                     %                         J_D..dbdphidphidphi = permute(sum(bsxfun(@times,J_D..dphidphidphidphi,permute(phi.db,[1,3,4,5,2])),1),[5,2,3,4,1]);
%                     %                         J_D..dbdbdphidphi = permute(sum(bsxfun(@times,J_D..dbdphidphidphi,permute(phi.db,[3,1,4,5,2])),2),[1,5,3,4,2]) ...
%                     %                             + permute(sum(bsxfun(@times,J_D.dphidphidphi,permute(phi.dbdb,[1,4,5,2,3])),1),[4,5,2,3,1]);
%                     %                         J_D..dbdbdbdphi = permute(sum(bsxfun(@times,J_D..dbdbdphidphi,permute(phi.db,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
%                     %                             + permute(sum(bsxfun(@times,J_D.dbdphidphi,permute(phi.dbdb,[4,1,5,2,3])),2),[1,4,5,3,2]);
%                     %                         J_D..dbdbdbdb = permute(sum(bsxfun(@times,J_D..dbdbdbdphi,permute(phi.db,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
%                     %                             + permute(sum(bsxfun(@times,J_D.dbdbdphi,permute(phi.dbdb,[4,5,1,2,3])),3),[1,2,4,5,3]);
%                     %
%                     %                         temp = squeeze(sum(bsxfun(@times,bsxfun(@times,J_T.dTdT,permute(T.dphidphi,[4,1,2,3])),permute(T.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_T.dTdR,permute(T.dphidphi,[4,1,2,3])),permute(R.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_T.dRdR,permute(R.dphidphi,[4,1,2,3])),permute(R.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_T.dTdSigma,permute(T.dphidphi,[4,1,2,3])),permute(Sigma_time.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_T.dRdSigma,permute(R.dphidphi,[4,1,2,3])),permute(Sigma_time.dphidphi,[4,1,5,6,2,3])) ...
%                     %                             + bsxfun(@times,bsxfun(@times,J_T.dSigmadSigma,permute(Sigma_time.dphidphi,[4,1,2,3])),permute(Sigma_time.dphidphi,[4,1,5,6,2,3])),2));
%                     %
%                     %                         J_T..dphidphidphidphi = ...
%                     %                             permute(temp,[1,3,2,4]) + permute(temp,[1,3,4,2]) + permute(temp,[3,4,1,2]);
%                     %
%                     %                         J_T..dphidTdTdT = bsxfun(@times,J_T..dTdTdTdT,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdTdTdR,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdTdTdSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidTdTdR = bsxfun(@times,J_T..dTdTdTdR,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdTdRdR,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdTdRdSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidTdRdR = bsxfun(@times,J_T..dTdTdRdR,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdRdRdR,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdRdRdSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidRdRdR = bsxfun(@times,J_T..dTdRdRdR,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dRdRdRdR,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dRdRdRdSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidTdTdSigma = bsxfun(@times,J_T..dTdTdTdSigma,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdTdRdSigma,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdTdSigmadSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidTdRdSigma = bsxfun(@times,J_T..dTdTdRdSigma,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdRdRdSigma,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdRdSigmadSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidRdRdSigma = bsxfun(@times,J_T..dTdRdRdSigma,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dRdRdRdSigma,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dRdRdSigmadSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidTdSigmadSigma = bsxfun(@times,J_T..dTdTdSigmadSigma,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdRdSigmadSigma,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdSigmadSigmadSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidRdSigmadSigma = bsxfun(@times,J_T..dTdRdSigmadSigma,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dRdRdSigmadSigma,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dTdSigmadSigmadSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidSigmadSigmadSigma = bsxfun(@times,J_T..dTdSigmadSigmadSigma,permute(T.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dRdSigmadSigmadSigma,permute(R.dphi,[3,1,2])) ...
%                     %                             + bsxfun(@times,J_T..dSigmadSigmadSigmadSigma,permute(Sigma_time.dphi,[3,1,2]));
%                     %                         J_T..dphidphidTdT = bsxfun(@times,J_T..dphidTdTdT,permute(T.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidTdTdR,permute(R.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidTdTdSigma,permute(Sigma_time.dphi,[3,1,4,2]));
%                     %                         J_T..dphidphidTdR = bsxfun(@times,J_T..dphidTdTdR,permute(T.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidTdRdR,permute(R.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidTdRdSigma,permute(Sigma_time.dphi,[3,1,4,2]));
%                     %                         J_T..dphidphidRdR = bsxfun(@times,J_T..dphidTdRdR,permute(T.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidRdRdR,permute(R.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidRdRdSigma,permute(Sigma_time.dphi,[3,1,4,2]));
%                     %                         J_T..dphidphidTdSigma = bsxfun(@times,J_T..dphidTdTdSigma,permute(T.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidTdRdSigma,permute(R.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidTdSigmadSigma,permute(Sigma_time.dphi,[3,1,4,2]));
%                     %                         J_T..dphidphidRdSigma = bsxfun(@times,J_T..dphidTdRdSigma,permute(T.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidRdRdSigma,permute(R.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidRdSigmadSigma,permute(Sigma_time.dphi,[3,1,4,2]));
%                     %                         J_T..dphidphidSigmadSigma = bsxfun(@times,J_T..dphidTdSigmadSigma,permute(T.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidRdSigmadSigma,permute(R.dphi,[3,1,4,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidSigmadSigmadSigma,permute(Sigma_time.dphi,[3,1,4,2]));
%                     %                         J_T..dphidphidphidT = bsxfun(@times,J_T..dphidphidTdT,permute(T.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidTdR,permute(R.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidTdSigma,permute(Sigma_time.dphi,[3,1,4,5,2]));
%                     %                         J_T..dphidphidphidR = bsxfun(@times,J_T..dphidphidTdR,permute(T.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidRdR,permute(R.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidRdSigma,permute(Sigma_time.dphi,[3,1,4,5,2]));
%                     %                         J_T..dphidphidphidSigma = bsxfun(@times,J_T..dphidphidTdSigma,permute(T.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidRdSigma,permute(R.dphi,[3,1,4,5,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidSigmadSigma,permute(Sigma_time.dphi,[3,1,4,5,2]));
%                     %                         J_T..dphidphidphidphi = J_T..dphidphidphidphi + squeeze(sum(bsxfun(@times,J_T..dphidphidphidT,permute(T.dphi,[3,1,4,5,6,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidphidR,permute(R.dphi,[3,1,4,5,6,2])) ...
%                     %                             + bsxfun(@times,J_T..dphidphidphidSigma,permute(Sigma_time.dphi,[3,1,4,5,6,2])),2));
%                     %
%                     %                         J_T..dbdphidphidphi = permute(sum(bsxfun(@times,J_T..dphidphidphidphi,permute(phi.db,[1,3,4,5,2])),1),[5,2,3,4,1]);
%                     %                         J_T..dbdbdphidphi = permute(sum(bsxfun(@times,J_T..dbdphidphidphi,permute(phi.db,[3,1,4,5,2])),2),[1,5,3,4,2]) ...
%                     %                             + permute(sum(bsxfun(@times,J_T.dphidphidphi,permute(phi.dbdb,[1,4,5,2,3])),1),[4,5,2,3,1]);
%                     %                         J_T..dbdbdbdphi = permute(sum(bsxfun(@times,J_T..dbdbdphidphi,permute(phi.db,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
%                     %                             + permute(sum(bsxfun(@times,J_T.dbdphidphi,permute(phi.dbdb,[4,1,5,2,3])),2),[1,4,5,3,2]);
%                     %                         J_T..dbdbdbdb = permute(sum(bsxfun(@times,J_T..dbdbdbdphi,permute(phi.db,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
%                     %                             + permute(sum(bsxfun(@times,J_T.dbdbdphi,permute(phi.dbdb,[4,5,1,2,3])),3),[1,2,4,5,3]);
%                     %
%                     %                         FIM.dbdb = J_D..dbdbdbdb + J_T..dbdbdbdb + squeeze(J_b..dbdbdbdb);
%                     %
%                     %                         %% FIM.dbdbeta
%                     %
%                     %                         J_D..dbdbdbdbeta = permute(sum(bsxfun(@times,J_D..dbdbdbdphi,permute(phi.dbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
%                     %                             + permute(sum(bsxfun(@times,J_D.dbdbdphi,permute(phi.dbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
%                     %
%                     %                         J_T..dbdbdbdbeta = permute(sum(bsxfun(@times,J_T..dbdbdbdphi,permute(phi.dbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
%                     %                             + permute(sum(bsxfun(@times,J_T.dbdbdphi,permute(phi.dbdbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
%                     %
%                     %                         FIM.dbdbeta = J_D..dbdbdbdbeta + J_T..dbdbdbdbeta;
%                     %
%                     %                         %% FIM.dbetadbeta
%                     %
%                     %                         J_D..dbdbdbetadphi = permute(sum(bsxfun(@times,J_D..dbdbdphidphi,permute(phi.dbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
%                     %                             + permute(sum(bsxfun(@times,J_D.dbdphidphi,permute(phi.dbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
%                     %                         J_D..dbdbdbetadbeta = permute(sum(bsxfun(@times,J_D..dbdbdbetadphi,permute(phi.dbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
%                     %                             + permute(sum(bsxfun(@times,J_D.dbdbdphi,permute(phi.dbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
%                     %
%                     %                         J_T..dbdbdbetadphi = permute(sum(bsxfun(@times,J_T..dbdbdphidphi,permute(phi.dbeta,[3,4,1,5,2])),3),[1,2,5,4,3]) ...
%                     %                             + permute(sum(bsxfun(@times,J_T.dbdphidphi,permute(phi.dbdbeta,[4,1,5,2,3])),2),[1,4,5,3,2]);
%                     %                         J_T..dbdbdbetadbeta = permute(sum(bsxfun(@times,J_T..dbdbdbetadphi,permute(phi.dbeta,[3,4,5,1,2])),4),[1,2,3,5,4]) ...
%                     %                             + permute(sum(bsxfun(@times,J_T.dbdbdphi,permute(phi.dbetadbeta,[4,5,1,2,3])),3),[1,2,4,5,3]);
%                     %
%                     %                         FIM.dbetadbeta = J_D..dbdbdbetadbeta + J_T..dbdbdbetadbeta;
%                     %
%                     %                         %% FIM.dbddelta
%                     %
%                     %                         FIM.dbddelta = J_b..dbdbdbddelta;
%                     %
%                     %                         %% FIM.ddeltaddelta
%                     %
%                     %                         FIM.ddeltaddelta = J_b..dbdbddeltaddelta;
%                     %
%                     %                         %% FIM.dbetaddelta
%                     %                         FIM.dbetaddelta = zeros(length(b),length(b),length(beta),length(b));
%                     
%                 end
            end
        end
    end
end

if nargout >=1
    varargout{1} = J.val;
end
if nargout >= 2
    varargout{2} = J.db;
end
if nargout >= 3
    varargout{3} = FIM.val;
end
if nargout >= 4
    varargout{4} = J;
    varargout{5} = FIM;
    varargout{6} = Sim;
end

end

function J = objective_phi_J_D(model,data,phi,s,i,options,nderiv)
    [J,~] = objective_phi(model,data,phi,s,i,options,nderiv);
end


function J = objective_phi_J_T(model,data,phi,s,i,options,nderiv)
    [~,J] = objective_phi(model,data,phi,s,i,options,nderiv);
end




