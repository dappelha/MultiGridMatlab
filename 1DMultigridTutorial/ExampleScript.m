%% 1-D Multigrid Code: Understanding the process 
%% Author: David Appelhans, Dec 22, 2011

% This code was written to better help you understand the multigrid process
% in the simple setting of the 1-D laplace equation with zeros boundary
% conditions,
%                  -u''(x) = f(x) in [a,b]      Differential Equation
%                   u(a) = 0, u(b) = 0.         Boundary Conditions
%
% You can examine the process of the multigrid solve by picking a u(x) that
% satisfies the boundary conditions and setting f(x) as the second
% derivative of the u(x) you have picked. Now you know what the exact
% solution u(x) is and you can see how close the numerical solution comes
% to this analytic solution.

% This code will be slow because we will generate plots at each step of the
% process.
%% Close all figures and clear all variables:
close all
clear all
%% Define differential equation parameters

% Specify the domain [a,b]:
a = 0; b = 1.5;

% Define the RHS function and the known analytic solution:
k = 2; L = b-a;
analyticV = @(x) sin(k*pi*x/L)';
f = @(x) (k*pi/L)^2*sin(k*pi*x/L)';

%% Define discretization parameters
% For this example I will be using the second order central difference
% formula: -u''(x_i) = (-u(x_{i+1}) + 2u(x_i) - u(x_{i-1}))/h^2 + O(h^2).
%
% At the end of the day this will lead to a linear algebra system Av=f
% where A is a matrix, f is a vector containing the values of the f(x)
% evaluated at a discrete grid of points and v is a vector of the
% approximate values of u(x) evaluated at the grid points x_i.

% The error from the discretization is O(h^2) meaning that our vector of
% function values v_i and the exact solution evaluated on the grid, u(x_i)
% may differ by about h^2. Thus there is no point in solving the
% linear algebra problem Av=f beyond O(h^2) since the linear algebra
% solution will only be accurate to O(h^2).

% Since the boundary conditions give us the values of u(x_i) on the
% boundaries, we only need to solve for the unknown function values in the
% interior of the grid

% Number of points to use in the numerical solution:

% Specify the number of interior points to use in the numerical solution
p = 4; % 9 gives 512 points
nx = 2^p -1 ; %(for multigrid this should be odd)
dx = (b-a)/(nx+1); % mesh spacing h=dx.

%% Multigrid Paramters
% Here you can change the number of grid levels, the type of multigrid
% cycle, etc.

% vis is a flag to visualize each step of the MuCycle. 
% vis = 1 for plots of the error approximations and corrections
% vis = 2 for plots of the residual on each grid level.
global vis
vis = 1;

% Specify Mu, Mu=1 vcycle, Mu = 2 wcycle
global Mu
Mu = 1;
% Specify the number of pre and post relaxations i.e. V(nu1,nu2).
nu1 = 2;
nu2 = 2;
% number of Mu cycles total to run:
nMucycles = 1;

% total number of grids:
global ngrids
ngrids = p-1;
%% 

[v] = VisualMultigrid1D(a,b,nx,nMucycles,nu1,nu2,analyticV,f);

%% Test weighted Jacobi:
% clear all
% close all
% % Specify the domain [a,b]:
% a = 0; b = 1.5;
% 
% % Define the RHS function and the known analytic solution:
% k = 0; L = b-a;
% exactV = @(x) sin(k*pi*x/L)';
% f = @(x) (k*pi/L)^2*sin(k*pi*x/L)';
% 
% p = 7; % 9 gives 512 points
% nx = 2^p -1 ; %(for multigrid this should be odd)
% dx = (b-a)/(nx+1); % mesh spacing h=dx.
% 
% xgrid = a+dx:dx:b-dx;
% ffunc = f(xgrid);
% ex = ones(nx,1); %array of ones
% A = spdiags([-ex 2*ex -ex], -1:1, nx, nx)/dx^2;
% vinitial = rand(nx,1);
% vfinal = WeightedJacobi(vinitial,ffunc,A,10000,2/3);
% figure
% plot(xgrid,vinitial,xgrid,vfinal,'--.')



