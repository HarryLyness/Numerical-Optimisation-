function [x,n] = bfgs_ex(f,df,df2,B,xn,tol)
% Generalized Steepest descent method (Alg. 4.3) with exact line
% search
%
% Inputs:
%
%   function f, gradient df, hessian df2
%   intitial guess xn  
%   tolerance for stopping criterion tol
%   initial B_0 matrix
%
% Outputs:
%
%   approximate solution x 
%   iterations n
%
% sets maximum iteration stopping criterion 
maxit=10000;
% computes the inverse of B_0
Binv = inv(B);

% pre-computes the correct dimension idenity matrix
I = eye(length(xn));

% start main loop, first stopping condition max iteration 
for n = 0:maxit
    x(:,n+1) = xn; 
    % Evaluate the Gradient at the current point x  
    gx = df(xn);
    % For the stopping criterion check whether the norm of the gradient is 
    % below the requested tolerance
    if norm(gx) <= tol   
        break
    end  
    % compute the decent direction at iteration n 
    sn = (-1)*Binv*gx;
    % Evaluate the hessian at the current point x
    df2x = df2(xn);
    % compute the alpha using exact linesearch method 
    alpha = exact_linesearch(gx,df2x,sn);   
    % update the solution 
    xn  = xn + alpha*sn;
    % update the inverse of B using the formula defined in (3) on
    % assignment 1
    dn = xn-x(:,n+1);
    yn = df(xn)-gx;
    rhon = 1/(transpose(yn)*dn);
    Binv = (I-rhon*dn*transpose(yn))*Binv*(I-rhon*yn*transpose(dn))+rhon*(dn*transpose(dn));
end
end