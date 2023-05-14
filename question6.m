%% question 6 function handles (RUN THIS)
clear
% Use function handles/ inline to define f, its gradient and Hessian
% gradient of f is a column vector df
% Hessian of f is a matrix d2f
f   = @(x) (x(1))^2 -x(1)*x(2)+5*(x(2))^2-2*x(1)+x(2);
df  = @(x) [2*x(1)-x(2)-2; (-1)*x(1)+10*x(2)+1]; 
d2f = @(x) [2, -1; -1, 10]; 
%% BFGS with exact line search with a range of initial x0 (RUN THIS)
% starting values not changed within each experiment
tol = 1.0e-5;
xex = [1;0];
B = eye(length(xex));
%% Function Random Starting values [-1.1,1.1] test: 
% Run 'R=10' many random tests, and define max and min values for range of 
% starting value x0
R = 10; 
min =-1.1; 
max = 1.1; 
% Do not need extra convergence infomation, so Info = 0 (more detail in
% dedicated functions e.g. 'R_X0_bfgs_ex.m')
Info = 0;
% computes the GSDM with Exact Linesearch for 10 random starting values in
% appropriate range and displays results in a table. 
disp(['GSDM with Exact Linesearch for ', num2str(R),' random starting ' ...
    'values in range ', num2str(min),' to ',num2str(max)])
Rtable = R_X0_bfgs_ex(f,df,d2f,B,tol,xex,min,max,R,Info);
disp(Rtable)

%% Function Random Starting values [-2,2] test: 
% Run 'R=10' many random tests, and define max and min values for range of 
% starting value x0
R = 10; 
min =-2; 
max = 2; 
% Do not need extra convergence infomation, so Info = 0 (more detail in
% dedicated functions e.g. 'R_X0_bfgs_ex.m')
Info = 0;
% computes the GSDM with Exact Linesearch for 10 random starting values in
% appropriate range and displays results in a table. 
disp(['GSDM with Exact Linesearch for ', num2str(R),' random starting ' ...
    'values in range ', num2str(min),' to ',num2str(max)])
Rtable = R_X0_bfgs_ex(f,df,d2f,B,tol,xex,min,max,R,Info);
disp(Rtable)
%% Function Random Starting values [-10,10] test: 
% Run 'R=10' many random tests, and define max and min values for range of 
% starting value x0
R = 10; 
min =-10; 
max = 10; 
% Do not need extra convergence infomation, so Info = 0 (more detail in
% dedicated functions e.g. 'R_X0_bfgs_ex.m')
Info = 0;
% computes the GSDM with Exact Linesearch for 10 random starting values in
% appropriate range and displays results in a table. 
disp(['GSDM with Exact Linesearch for ', num2str(R),' random starting ' ...
    'values in range ', num2str(min),' to ',num2str(max)])
Rtable = R_X0_bfgs_ex(f,df,d2f,B,tol,xex,min,max,R,Info);
disp(Rtable)
%% Description of following non-random starting value tests
% Here are some examples for non-random starting values. The first test
% trajectory is included in the written solutions pdf file. 
% included solution trajectories for these because they are pretty :) 
%% TEST1 (trajectory in written solutions pdf file)
% Define starting value x0 and compute convergence tests using dedicated
% function below for bfgs with exact line search 
x0 = [-1;1];
[xsol,nsol] = bfgs_ex_convergence(f,df,d2f,B,x0,tol,xex);
%Solution path for Generalised Steepest Decent with Exact Linesearch
%Algorithm with f(x) = (1/2)x^TAx + b^Tx and starting value x_0 = [-1,1]^T
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%% TEST2 (x0 close to minimum)
% Define starting value x0 and compute convergence tests using dedicated
% function below for bfgs with exact line search 
x0 = [0.9;1.1];
[xsol,nsol] = bfgs_ex_convergence(f,df,d2f,B,x0,tol,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%% TEST3 (x0 further away from minimum)
% Define starting value x0 and compute convergence tests using dedicated
% function below for bfgs with exact line search 
x0 = [-2;-3];
[xsol,nsol] = bfgs_ex_convergence(f,df,d2f,B,x0,tol,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%% TEST4 (x0 far away from minimum)
% Define starting value x0 and compute convergence tests using dedicated
% function below for bfgs with exact line search 
x0 = [50;100];
[xsol,nsol] = bfgs_ex_convergence(f,df,d2f,B,x0,tol,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
