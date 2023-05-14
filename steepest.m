function [x,n] = steepest(f,df,xn,theta,tol,maxit)
% Steepest descent method (Alg. 4.2) with backtracking line search (Alg. 4.1)
%
% Inputs:
%
%   function f and gradient df
%   intitial guess xn
%   user defined theta for Armijo condition
%   tolerance for stopping criterion tol
%
% Output:
%   approximate solution x 
%   iterations n


% Iterate up to a maximum of maxit iterations
for n = 0:maxit

    x(:,n+1) = xn;

    % Evaluate the Gradient at the current point x    
    gx = df(xn);

    % For the stopping criterion check whether the norm of the gradient is 
    % below the requested tolerance
    if norm(gx) < tol   
        break
    end
    
    % Evaluate the objective function at the current point x    
    fx = f(xn);

    % Call backtracking line search with steepest descent direction    
    alpha = linesearch(xn,fx,gx,-gx,theta,f); 
    
    % Calculate the next iterate in steepest descent direction    
    xn  = xn - alpha*gx;             
end

end