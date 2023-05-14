function [x,n] = bfgs(f,df,B,xn,theta,tol)
% Generalized Steepest descent method (Alg. 4.3) with backtracking line
% search (Alg. 4.1)
%
% Inputs:
%
%   function f and gradient df
%   intitial guess xn
%   user defined theta for Armijo condition
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
    % Evaluate the objective function at the current point x 
    fx = f(xn);  
    % compute the decent direction at iteration n 
    sn = (-1)*Binv*gx;
    % compute the alpha using linesearch method 
    alpha = linesearch(xn,fx,gx,sn,theta,f);    
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