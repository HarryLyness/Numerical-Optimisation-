function [xsol,nsol] = newtonconvergence(df,d2f,x0,tol,xex)
% FUNCTION to compute rates of convergence, number of iterations, and
% final error for newtons method
%
% INPUTS: 
%         df - gradient of f
%         d2f - hessian of f
%         x0 - starting value of x
%         tol - ubser defined tolerance 
%         xex - exact solution
%         
% OUTPUTS:
%        xsol - computed solution x_n using newtons method  
%        nsol - number of iterations x_0,...,x_n 
%         
%         
% compute newtons method for starting value x0.
[xsol,nsol] = newton(df,d2f,x0,tol);
%display number of iterations
disp(['Number of iterations: ', num2str(nsol)]);
% make testing convergence tables 
disp(['Table for Newtons method for x0 =[',num2str(transpose(x0)),']^T'])
disp(['Table designed to show potential candidate convergence rates'])
[output,e] = makeconvergencetable(xex,xsol,nsol);
disp(output)
% displays final error
disp(['Final error: ', num2str(e(nsol+1))]);
end 