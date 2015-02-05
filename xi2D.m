function [D,invD,dD,dinvD,HD,HinvD] = xi2D(xi,type)

% Covariance matrix D
switch type
    case 'diag-matrix-logarithm'
        % Diagonal matrix logarithm parameterization
        m = length(xi);
        D = diag(exp(xi));
        invD = diag(exp(-xi));
        dD = zeros(m,m,m);
        dinvD = zeros(m,m,m);
        HD = zeros(m,m,m,m);
        HinvD = zeros(m,m,m,m);
        for i = 1:length(xi)
            dD(i,i,i) = exp(xi(i));
            dinvD(i,i,i) = -exp(-xi(i));
            HD(i,i,i,i) = exp(xi(i));
            HinvD(i,i,i,i) = exp(-xi(i));
        end
    case 'matrix-logarithm'
        % Matrix logarithm parameterization
        m = round((-1 + sqrt(1+8*length(xi)))/2);
        logD = zeros(m,m);
        dlogD = zeros(m,m,length(xi));
        d2logD = zeros(m,m,length(xi),length(xi));
        dD = zeros(m,m,length(xi));
        dinvD = zeros(m,m,length(xi));
        HD = zeros(m,m,length(xi),length(xi));
        HinvD = zeros(m,m,length(xi),length(xi));
        
        k = 1;
        for i = 1:m
            for j = 1:i
                logD(i,j) = xi(k);
                logD(j,i) = xi(k);
                dlogD(i,j,k) = 1;
                dlogD(j,i,k) = 1;
                k = k+1;
            end
        end
        
        for i = 1:length(xi)
            for j = 1:i
                M = expm([logD(:,:)   , dlogD(:,:,i) , dlogD(:,:,j) , d2logD(:,:,i,j) ;...
                          zeros(m)    , logD(:,:)    , zeros(m)     , dlogD(:,:,j) ;...
                          zeros(m)    , zeros(m)     , logD(:,:)    , dlogD(:,:,i) ;...
                          zeros(m)    , zeros(m)     , zeros(m)     , logD(:,:)]);
                invM = expm([-logD(:,:)   , -dlogD(:,:,i) , -dlogD(:,:,j) , -d2logD(:,:,i,j) ;...
                          zeros(m)    , -logD(:,:)    , zeros(m)     , -dlogD(:,:,j) ;...
                          zeros(m)    , zeros(m)     , -logD(:,:)    , -dlogD(:,:,i) ;...
                          zeros(m)    , zeros(m)     , zeros(m)     , -logD(:,:)]);
                HD(:,:,i,j) = M(1:m,3*m+1:4*m);
                HD(:,:,j,i) = HD(:,:,i,j)';
                HinvD(:,:,i,j) = invM(1:m,3*m+1:4*m);
                HinvD(:,:,j,i) = HinvD(:,:,i,j)';
            end
            dD(:,:,i) = M(1:m,m+1:2*m);
            dinvD(:,:,i) = invM(1:m,m+1:2*m);
        end
        D = M(1:m,1:m);
        invD = invM(1:m,1:m);
        
        
%         D = expm(logD);
%         invD = expm(-logD);
%         
%         [X,L] = eig(logD);
%         invX = inv(X);
%         V = zeros(m,m);
%         for i = 1:m
%             for j = 1:m
%                 if i == j
%                     V(i,i) = exp(L(i,i));
%                     invV(i,i) = exp(-L(i,i));
%                 else
%                     V(i,j) = (exp(L(i,i))-exp(L(j,j)))/(L(i,i)-L(j,j));
%                     invV(i,j) = (exp(-L(i,i))-exp(-L(j,j)))/(-L(i,i)+L(j,j));
%                 end
%             end
%         end
%         
%         dD = zeros(m,m,length(xi));
%         dinvD = zeros(m,m,length(xi));
%         for i = 1:length(xi)
%             G = invX*dlogD(:,:,i)*X;
%             dD(:,:,i) = X*(V.*G)*invX;
%             dinvD(:,:,i) = X*(-invV.*G)*invX;
%         end

%         eps = 1e-5; 
%         for k = 1:length(xi)
%             xi_k = xi;
%             xi_l = xi;
%             xi_k(k) = xi_k(k) + eps;
%             xi_l(k) = xi_l(k) - eps;
%             logD_k = [xi_k(1),xi_k(2),xi_k(4),xi_k(7);
%                       xi_k(2),xi_k(3),xi_k(5),xi_k(8);
%                       xi_k(4),xi_k(5),xi_k(6),xi_k(9);
%                       xi_k(7),xi_k(8),xi_k(9),xi_k(10)];
%             dDdxi_fd(:,:,k) = (expm(logD_k) - D)/eps;
%             for l = 1:length(xi)
%                 xi_kl = xi_k;
%                 xi_lk = xi_l;
%                 xi_kl(l) = xi_kl(l) + eps;
%                 xi_lk(l) = xi_lk(l) - eps;
%                 logD_kl = [xi_kl(1),xi_kl(2),xi_kl(4),xi_kl(7);
%                            xi_kl(2),xi_kl(3),xi_kl(5),xi_kl(8);
%                            xi_kl(4),xi_kl(5),xi_kl(6),xi_kl(9);
%                            xi_kl(7),xi_kl(8),xi_kl(9),xi_kl(10)];
%                 logD_lk = [xi_lk(1),xi_lk(2),xi_lk(4),xi_lk(7);
%                            xi_lk(2),xi_lk(3),xi_lk(5),xi_lk(8);
%                            xi_lk(4),xi_lk(5),xi_lk(6),xi_lk(9);
%                            xi_lk(7),xi_lk(8),xi_lk(9),xi_lk(10)];
%                        if(k==l)
%                            HDdxi_fd(:,:,k,l) = (expm(logD_kl) + expm(logD_lk) - 2*D)/(4*eps^2);
%                        else
%                            HDdxi_fd(:,:,k,l) = (expm(logD_kl) + expm(logD_lk) - 2*D)/eps;
%                        end
%             end
%             
%         end
%         
%         D
%         dD
%         dD - dDdxi_fd
%         HD - HDdxi_fd
        
end


