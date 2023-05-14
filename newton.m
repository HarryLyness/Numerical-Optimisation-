function [x,n] = newton(df,d2f,xn,tol)
% Newton method for finding gradient df = 0
%
% Inputs:
%
%   gradient df, Hessian d2f
%   intitial guess xn
%   tolerance for stopping criterion tol
%
% Output:
%   approximate solution x 
%   iterations n
%
% Iterate up to a maximum of 100 iterations

for n=0:100
    
    x(:,n+1) = xn;

    % Evaluate the Gradient at the current point x
    
    gx = df(xn);

    % For the stopping criterion check whether the norm of the gradient is 
    % below the requested tolerance
    if norm(gx) < tol
        break
    end
    
    % Evaluate the Hessian at the current point x and calculate the 
    % next iterate via Newton's Method    
    B = d2f(xn);
    xn = xn - B \ gx;
        
end

% Check that the critical point that was found is indeed a local minimum

if min(eig(B)) > 0 
    disp('Found a local minimum.');
else
    disp('Newton ended up in a critical point that is not a local minimum');
end

end