%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR COMPUTATION OF DERIVATIVES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdz = chainrule_dxdy_dydz(dxdy,dydz)
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
            dxdz = permute(sum(bsxfun(@times,dxdy,permute(dydz,[d2+(1:(d1-1)),1,2:d2])),d1),[1:(d1-1),(d1+1):(d1+d2-1),d1]);
        end
    else
        dxdz =     permute(sum(bsxfun(@times,dxdy,permute(dydz,[              1,2:d2])),d1),[1:(d1-1),(d1+1):(d1+d2-1),d1]);
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

function ddxdzdz = chainrule_ddxdydy_dydz(ddxdydy,dydz)
d1 = ndims(ddxdydy);
d2 = ndims(dydz);
%      dx        dy      dy
% [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d2-1]
%               <.,.>
%   ********     dy     ***           dz       **************
% [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d2-1]
%   ******************   dy      **********         dz
%                       <.,.>
% [ 1 : d1 - 1        , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d2-1]
if(d1>2)
    ddxdzdy =             sum(bsxfun(@times,ddxdydy,permute(dydz,[d2+(1:(d1-2)),1,d2+d1-1,2:d2             ]     )),d1-1);
else
    if(size(dydz,1)==1)
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[                d2+d1-1,1:d2             ]     )),d1);
    else
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[              1,d2+d1-1,2:d2             ]     )),d1-1);
    end
end
ddxdzdz =         squeeze(sum(bsxfun(@times,ddxdzdy,permute(dydz,[d2+(1:(d1-1))  ,1      ,d2+(d1:(d1+d2-1)),2:d2])),d1));

end

function ddxdzdv = chainrule_ddxdydy_dydz_dydv(ddxdydy,dydz,dydv)
d1 = ndims(ddxdydy);
d2 = ndims(dydz);
d3 = ndims(dydv);
%      dx        dy      dy
% [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d3-1]
%               <.,.>
%   ********     dy     ***           dz       **************
% [ 1 : d1-2 , d1 - 1 , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d3-1]
%   ******************   dy      **********         dv
%                       <.,.>
% [ 1 : d1 - 1        , d1     , d1 + 1 : d2-1 , d1 + d2 + 1:d3-1]
if(d1>2)
    ddxdzdy =             sum(bsxfun(@times,ddxdydy,permute(dydz,[d2+(1:(d1-2)),1,d2+d1-1,2:d2             ]     )),d1-1);
else
    if(size(dydz,1)==1)
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[                d2+d1-1,1:d2             ]     )),d1);
    else
        ddxdzdy =         sum(bsxfun(@times,ddxdydy,permute(dydz,[              1,d2+d1-1,2:d2             ]     )),d1-1);
    end
end
ddxdzdv =         squeeze(sum(bsxfun(@times,ddxdzdy,permute(dydv,[d3+(1:(d1-1))  ,1      ,d3+(d1:(d1+d2-1)),2:d3])),d1));

end