%% question 5 (RUN THIS)
clear
% Use function handles/ inline to define f, its gradient and Hessian
% gradient of f is a column vector df
% Hessian of f is a matrix d2f
f   = @(x) (1-x(1))^2 + 100 * (x(2) -x(1)^2)^2;
df  = @(x) [2*(x(1)-1) + 400*x(1)*(x(1)^2 - x(2)); 200*(x(2) - x(1)^2)];
d2f = @(x) [2-400*x(2)+1200*x(1)^2, -400*x(1); -400*x(1), 200];           

%% generalised steepest descent method with wlinesearch for x0 =[-1.2;1]^T
% Initial guess
x0  = [-1.2;1];

% Necessary inputs
tol = 1.0e-5;
theta_sd= 0.1;
theta_c = 0.9;
B = eye(length(x0));

% exact solution computed analytically
xex = [1;1];

% compute x_n solution and store in two arrays 
[xsol,nsol] = bfgsw(f,df,B,x0,theta_sd,tol,theta_c);

% Visualise the solution path in the contour plot
%BFGS with Wolfe Linesearch Algorithm for Rosenbrock 
% function with starting value x0 = [-1.2,1]^T
visual(f,xsol,x0,xex);

%display number of iterations
disp(['Number of iterations: ', num2str(nsol)]);

% make testing convergence tables 
disp(['Table for generalised steepest descent method with wlinesearch for x0 =[-1.2;1]^T'])
disp(['Table designed to show potential candidate convergence rates'])
[output,e] = makeconvergencetable(xex,xsol,nsol);
disp(output)

% displays final error
% Check accuracy of the solution
disp(['Final error: ', num2str(e(nsol+1))]);

%% generalised steepest descent method with wlinesearch for x0 =[-1.2;1.5]^T
% Initial guess
x0  = [-1.2;1.5];
% Necessary inputs
tol = 1.0e-5;
theta_sd= 0.1;
theta_c = 0.9;
B = eye(length(x0));
% exact solution computed analytically
xex = [1;1];
% compute x_n solution and store in two arrays 
[xsol,nsol] = bfgsw(f,df,B,x0,theta_sd,tol,theta_c);
%BFGS with Wolfe Linesearch Algorithm for Rosenbrock 
%function with starting value x0 = [-1.2,1.5]^T
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%display number of iterations
disp(['Number of iterations: ', num2str(nsol)]);
% make testing convergence tables 
disp(['Table for generalised steepest descent method with wlinesearch for x0 =[-1.2;1.5]^T'])
disp(['Table designed to show potential candidate convergence rates'])
[output,e] = makeconvergencetable(xex,xsol,nsol);
disp(output)
% displays final error
% Check accuracy of the solution
disp(['Final error: ', num2str(e(nsol+1))]);

%% Creating tables and plots for NEWTON with x0 = [-1.2;1]^T
% initialise starting values for Newton
x0  = [-1.2;1];
tol = 1.0e-5;
% Set exact solution
xex = [1;1];
%compute rates of convergence, number of iterations, and
% final error for newtons method
[xsol,nsol] = newtonconvergence(df,d2f,x0,tol,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%% Creating tables and plots for NEWTON with x0  = [-1.2;1.5]^T
% initialise starting values for Newton
x0  = [-1.2;1.5];
tol = 1.0e-5;
% Set exact solution
xex = [1;1];
%compute rates of convergence, number of iterations, and
% final error for newtons method
[xsol,nsol] = newtonconvergence(df,d2f,x0,tol,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%% Creating tables and plots for Steepest decent with linesearch and x0 = [-1.2;1]^T
% initialise starting values for Steepest decent
x0  = [-1.2;1];
tol = 1.0e-5;
theta=0.1;
maxit = 10000;
% Set exact solution
xex = [1;1];
%compute rates of convergence, number of iterations, and
% final error for steepest decent method with linesearch
[xsol,nsol] = steepestconvergence(f,df,x0,theta,tol,maxit,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);
%% Creating tables and plots for Steepest decent with linesearch and x0  = [-1.2;1.5]^T
% initialise starting values for Steepest decent
x0  = [-1.2;1.5];
tol = 1.0e-5;
theta=0.1;
maxit = 10000;
% Set exact solution
xex = [1;1];
%compute rates of convergence, number of iterations, and
% final error for steepest decent method with linesearch
[xsol,nsol] = steepestconvergence(f,df,x0,theta,tol,maxit,xex);
% Visualise the solution path in the contour plot
visual(f,xsol,x0,xex);