%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTION FOR DERIVATIVE OF OPTIMA W.R.T HYPERPARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [bhat,dbhatdbeta,dbhatddelta,ddbhatdbetadbeta,ddbhatdbetaddelta,ddbhatddeltaddelta] = bhat_SCTL_si(bhat,G,pdGpdbeta,pdGpddelta);
function varargout = bhat_SCTL_si(bhat,G,ddJdbdbeta,ddJdbddelta,~,~,~,dGdb,pdGpdbeta,pdGpddelta,dddJdbdbetadbeta,dddJdbddeltaddelta,dddJdbdbetaddelta)

invG = pinv(G);

if nargout>=2
    
    dbhatdbeta = -invG*squeeze(ddJdbdbeta);
    dbhatddelta = -invG*squeeze(ddJdbddelta);
    
    if nargout >= 4
        
        ddbhatdbetadbeta = zeros(size(dddJdbdbetadbeta));
        ddbhatdbetaddelta = zeros(size(dddJdbdbetaddelta));
        ddbhatddeltaddelta = zeros(size(dddJdbddeltaddelta));
        
        for j=1:length(bhat)
            for k = 1:length(bhat)
                if length(bhat) == 1
                    ddbhatdbetadbeta = squeeze(ddbhatdbetadbeta) - invG*(...
                        squeeze(dddJdbdbetadbeta) ...
                        + transpose(dbhatdbeta)*transpose(squeeze(pdGpdbeta(k,:,:))) ...
                        + transpose(transpose(dbhatdbeta)*transpose(squeeze(pdGpdbeta(k,:,:)))) ...
                        + transpose(dbhatdbeta)*squeeze(dGdb(k,:,:))*dbhatdbeta);
                    ddbhatdbetadbeta = permute(ddbhatdbetadbeta,[3,1,2]);
                    
                    ddbhatdbetaddelta = transpose(squeeze(ddbhatdbetaddelta)) - invG*(...
                        transpose(squeeze(dddJdbdbetaddelta)) ...
                        + transpose(dbhatdbeta)*squeeze(pdGpddelta(k,:,:)) ...
                        + transpose(dbhatddelta)*squeeze(pdGpdbeta(k,:,:)) ...
                        + transpose(dbhatdbeta)*squeeze(dGdb(k,:,:))*dbhatddelta);
                    ddbhatdbetaddelta = permute(ddbhatdbetaddelta,[3,1,2]);
                    
                    ddbhatddeltaddelta = squeeze(ddbhatddeltaddelta(j,:,:)) -invG(j,k)*(...
                        squeeze(dddJdbddeltaddelta(k,:,:)) ...
                        + transpose(dbhatddelta)*squeeze(pdGpddelta(k,:,:)) ...
                        + transpose(transpose(dbhatddelta)*squeeze(pdGpddelta(k,:,:))) ...
                        + transpose(dbhatddelta)*squeeze(dGdb(k,:,:))*dbhatddelta);
                    ddbhatddeltaddelta = permute(ddbhatddeltaddelta,[3,1,2]);
                else
                    ddbhatdbetadbeta(j,:,:) = squeeze(ddbhatdbetadbeta(j,:,:)) - invG(j,k)*(...
                        squeeze(dddJdbdbetadbeta(k,:,:)) ...
                        + transpose(dbhatdbeta)*squeeze(pdGpdbeta(k,:,:)) ...
                        + transpose(transpose(dbhatdbeta)*squeeze(pdGpdbeta(k,:,:))) ...
                        + transpose(dbhatdbeta)*squeeze(dGdb(k,:,:))*dbhatdbeta);
                    
                    ddbhatdbetaddelta(j,:,:) = squeeze(ddbhatdbetaddelta(j,:,:)) - invG(j,k)*(...
                        squeeze(dddJdbdbetaddelta(k,:,:)) ...
                        + transpose(dbhatdbeta)*squeeze(pdGpddelta(k,:,:)) ...
                        + transpose(transpose(dbhatddelta)*squeeze(pdGpdbeta(k,:,:))) ...
                        + transpose(dbhatdbeta)*squeeze(dGdb(k,:,:))*dbhatddelta);
                    
                    ddbhatddeltaddelta(j,:,:) = squeeze(ddbhatddeltaddelta(j,:,:)) -invG(j,k)*(...
                        squeeze(dddJdbddeltaddelta(k,:,:)) ...
                        + transpose(dbhatddelta)*squeeze(pdGpddelta(k,:,:)) ...
                        + transpose(transpose(dbhatddelta)*squeeze(pdGpddelta(k,:,:))) ...
                        + transpose(dbhatddelta)*squeeze(dGdb(k,:,:))*dbhatddelta);
                end
            end
        end
    end
end



if nargout <= 1
    varargout{1} =  bhat;
elseif nargout <= 3
    varargout{1} =  bhat;
    varargout{2} =  dbhatdbeta;
    varargout{3} =  dbhatddelta;
elseif nargout <= 6
    varargout{1} =  bhat;
    varargout{2} =  dbhatdbeta;
    varargout{3} =  dbhatddelta;
    varargout{4} =  ddbhatdbetadbeta;
    varargout{5} =  ddbhatdbetaddelta;
    varargout{6} =  ddbhatddeltaddelta;
end
end