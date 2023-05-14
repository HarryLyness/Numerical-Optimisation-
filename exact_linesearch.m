function alpha = exact_linesearch(dfx,df2x,sn)
% function to compute alpha using exact linesearch
%
% INPUTS: 
%
%   dfx - gradient at current point x
%   df2x - hessian at current point x
%   sn - decent direction at current iteration n
%
% OUTPUT: 
%
%   alpha 
%
% compute alpha for exact line search using formula in problem sheet 3
alpha = ((-1)*transpose(dfx)*sn)/(transpose(sn)*df2x*sn);
end