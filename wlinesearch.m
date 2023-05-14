function alpha = wlinesearch(f,df,sn,xn,theta_sd,theta_c)
% function: implementation of wlinesearch algorithm 6.1 in the lecture notes
%
% inputs:    xn - current iterate
%           f - f(x)
%           df - gradient f(x)
%           sn - search direction 
%           theta_sd - user defined threshold from Armijo condition
%           theta_c - user defined curvature condition...
%
% outputs:   alpha
%
% initialise \overline{alpha}, \underline{alpha} and alpha
alpha =1; alpha_upper = 0; alpha_lower = 0;
% complete procedure while either the curvature condition or Armijo
% condition fails 
while (f(xn+alpha*sn) > f(xn) + theta_sd*alpha*dot(df(xn),sn)) || (dot(df(xn+alpha*sn),sn)<theta_c*transpose(df(xn))*sn)
    % does the Armijo condition fail?
    if (f(xn+alpha*sn) > f(xn) + theta_sd*alpha*dot(df(xn),sn))
        % re-initialise \overline{alpha}, \underline{alpha} and alpha
        % reduce alpha
        alpha_upper = alpha; alpha = 1/2*(alpha_upper+alpha_lower);
    % does the curvature condition fail?
    elseif (dot(df(xn+alpha*sn),sn)<theta_c*transpose(df(xn))*sn)
        % re-initialise \underline{alpha}
        % increase alpha
        alpha_lower = alpha;
        % alpha is a real number so adjust for MATLAB floating point errors
        if alpha_upper<10^(-15) && -10^(-15)<alpha_upper
            alpha = 2*alpha_lower;
        else
            alpha = 1/2*(alpha_lower+alpha_upper);
        end 
    end
end 
