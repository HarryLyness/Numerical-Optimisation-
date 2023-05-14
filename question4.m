%% question 4 (RUN THIS)
clear
% initiate starting x here for testing multiple starting values as required
% im lazy!
x0 = [-1;1];
%x0  = [-1.2;1]; %x0 = [-1;0.5]; %x0 = [-50;1000];

% Use function handles/ inline to define f, its gradient and Hessian
% gradient of f is a column vector df
% Hessian of f is a matrix d2f
f   = @(x) (x(1))^2 -x(1)*x(2)+5*(x(2))^2-2*x(1)+x(2);
df  = @(x) [2*x(1)-x(2)-2; (-1)*x(1)+10*x(2)+1]; 
d2f = @(x) [2, -1; -1, 10];  

%% Creating tables and plots for BFGS
% initialise starting values for bfgs
%x0 = [-1;1];
tol = 1.0e-5;
theta=0.1;
B = eye(length(x0));
% compute x_n solution and store in two arrays 
[xsol,nsol] = bfgs(f,df,B,x0,theta,tol);
% Set exact solution
xex = [1;0];
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%display number of iterations
disp(['Number of iterations: ', num2str(nsol)]);
% make testing convergence tables 
disp(['Table for generalised steepest descent method with backtracking' ...
    ' line search x0 =[-1;1]^T'])
disp(['Table designed to show potential candidate convergence rates'])
[output,e] = makeconvergencetable(xex,xsol,nsol);
disp(output)
% displays final error
disp(['Final error: ', num2str(e(nsol+1))]);

%% Creating tables and plots for NEWTON
% initialise starting values for Newton
%x0 = [-1;1];
tol = 1.0e-5;
% Set exact solution
xex = [1;0];
%compute rates of convergence, number of iterations, and
% final error for newtons method
[xsol,nsol] = newtonconvergence(df,d2f,x0,tol,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);

%% Creating tables and plots for Steepest decent with linesearch
% initialise starting values for steepest decent
%x0 = [-1;1];
tol = 1.0e-5;
theta=0.1;
maxit = 10000;
% Set exact solution
xex = [1;0];
%compute rates of convergence, number of iterations, and
% final error for steepest decent method with linesearch
[xsol,nsol] = steepestconvergence(f,df,x0,theta,tol,maxit,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);