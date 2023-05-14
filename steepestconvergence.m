function [xsol,nsol] = steepestconvergence(f,df,x0,theta,tol,maxit,xex)
% FUNCTION to compute rates of convergence, number of iterations, and
% final error for steepest decent method with linesearch
%
% INPUTS: 
%         df - gradient of f
%         f - function f defined by user
%         x0 - starting value of x
%         tol - ubser defined tolerance 
%         xex - exact solution
%         maxit - maximum iteration stopping condition defined by user 
%         theta - user defined threshold from Armijo condition
%
% OUTPUTS:        
%        xsol - computed solution x_n using newtons method  
%        nsol - number of iterations x_0,...,x_n 
%         
%         
% compute steepest decent method for starting value x0.
% compute x_n solution and store in two arrays 
[xsol,nsol] = steepest(f,df,x0,theta,tol,maxit);
%display number of iterations
disp(['Number of iterations: ', num2str(nsol)]);
% make testing convergence tables 
disp(['Table for Steepest decent with backtracking linesearch for x0=[',num2str(transpose(x0)),']^T'])
disp(['Table designed to show potential candidate convergence rates'])
[output,e] = makeconvergencetable(xex,xsol,nsol);
disp(output)
% displays final error
disp(['Final error: ', num2str(e(nsol+1))]);
end
