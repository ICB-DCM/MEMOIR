%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR COMPUTATION OF DERIVATIVES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generalisation of product A*b for higher order tensors

function dxdz = chainrule(dxdy,dydz)
d1 = ndims(dxdy);
d2 = ndims(dydz);

if(~size(dxdy,d1)==size(dydz,1))
    error('dimensions must agree')
end
%      dx        dy
% [ 1 : d1 - 1 , d1     , d1 + 1:d2-1]
%               <.,.>
%   ********     dy          dz
% [ 1 : d1 - 1 , d1     , d1 + 1:d2-1]
if(numel(dxdy)>1)
    if(d1>1)
        if(and(d1==2,d2==2))
            if(size(dydz,1)==1)
                dxdz = bsxfun(@times,dxdy,permute(dydz,[3,1,2]));
            else
                dxdz = dxdy*dydz;
            end
        else
            dxdz = permute(nansum(bsxfun(@times,dxdy,permute(dydz,[d2+(1:(d1-1)),1,2:d2])),d1),[1:(d1-1),(d1+1):(d1+d2-1),d1]);
        end
    else
        dxdz =     permute(nansum(bsxfun(@times,dxdy,permute(dydz,[              1,2:d2])),d1),[1:(d1-1),(d1+1):(d1+d2-1),d1]);
    end
elseif(numel(dxdy)==1)
    dxdz = permute(dxdy*dydz,[d2+(1:(d1-1)),2:d2,1]);
else
    s1 = size(dxdy);
    s2 = size(dydz);
    dxdz = zeros([s1(1:(d1-1)),s2(2:d2)]);
end
dxdz = squeeze(dxdz);
end