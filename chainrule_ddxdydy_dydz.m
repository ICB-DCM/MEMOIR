%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS FOR COMPUTATION OF DERIVATIVES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generalisation of product b'*A*b for higher order tensors

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