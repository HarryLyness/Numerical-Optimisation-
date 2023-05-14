function Rtable = R_X0_bfgs_ex(f,df,d2f,B,tol,xex,min,max,R,Info)
% function to compute Generalized Steepest descent method (Alg. 4.3) with exact line
% search for random starting values x0 in range of [min,max], 
% with option to generate 
%  - convergence information: number of iterations, convergence table and final error
%    for each randomly chosen x0 
%
% Inputs:
%
%   function f, gradient df, hessian df2
%   tolerance for stopping criterion tol
%   initial B_0 matrix
%   exact solution for minimum 'xex'
%   Range of vector values for initial starting value can be set by using 
%       - min: minimum value for starting value 
%       - max: maximum value for starting value
%       Note that for best results, follow convention min < max (both real
%       valued)
%   R - Number of random starting values to be created and displayed
%   Info - select '1' if you want to see the convergence table and final
%       error for each starting value x0. If not then select '0', which is
%       recommended if R is large (Reduce output to user). 
%   
% Outputs:
%
%    Rtable - table with R many starting values x0 with number of
%    iterations to converge to solution 
%    
% length of vector x0
N = length(xex);
% indicator variable to display if iterations till solution is greater than
% N
flag = 0;
for n = 1:R
    % generates random starting vector value x0
    x0 = min + (max-min).*rand(N,1);
    % saves starting value for table
    x0R(:,n) = x0;
    % if info is 1, then computes extra convergence information for each
    % starting value and displays infomation to user
    if Info == 1
        [xsol,nsol] = bfgs_ex_convergence(f,df,d2f,B,x0,tol,xex);
    else 
        [xsol,nsol] = bfgs_ex(f,df,d2f,B,x0,tol);
    end
    %saves # of iterations until solution for table
    nsolR(:,n) = nsol; 
    % is the number of iterations until solution > N? (test)
    if nsol > N
        % change flag indicator variable if so
        flag = 1;
    end 
end 
% if any of the iterations were over N, then displays message for user at
% bottom of table, this is useful for large R
if flag == 1
    disp('A random starting value x0 converged in more than 2 iterations')
end 
% creating table described above with appropriate titles
Rtable = table(num2str(x0R'),nsolR');
Rtable = renamevars(Rtable,['Var1'],['Starting value x_0']);
Rtable = renamevars(Rtable,['Var2'],['Number of iterations to reach solution']);

