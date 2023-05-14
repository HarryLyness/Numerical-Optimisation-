function [xsol,nsol] = bfgs_ex_convergence(f,df,d2f,B,x0,tol,xex)
% function to compute Generalized Steepest descent method (Alg. 4.3) with exact line
% search convergence information: number of iterations, convergence table
% and final error 
%
% Inputs:
%
%   function f, gradient df, hessian df2
%   intitial guess x0  
%   tolerance for stopping criterion tol
%   initial B_0 matrix
%   exact solution for minimum 'xex'
%
% Outputs:
%
%    xsol - computed iterates x_n
%    sol - number of iterates x_0,...,x_n
%
% compute x_n solution and store in two arrays 
[xsol,nsol] = bfgs_ex(f,df,d2f,B,x0,tol);
%display number of iterations
disp(['Number of iterations: ', num2str(nsol)]);
% make testing convergence tables 
disp(['Table for BFGS method with exact line search for x0 =[', num2str(transpose(x0)) ,']^T'])
disp(['Table designed to show potential candidate convergence rates'])
[output,e] = makeconvergencetable(xex,xsol,nsol);
disp(output)
% displays final error
disp(['Final error: ', num2str(e(nsol+1))]);
end 