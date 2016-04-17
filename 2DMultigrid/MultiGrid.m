function [pv,resnorm,errornorm] = MultiGrid(ax,bx,ay,by,nxfine,nyfine,nMucycles,nu1,nu2,exactV,f)
%% Author: David Appelhans 2011
%% solves the finite difference discretization of square 2-D Poisson.
% Specify the number of grids (must be <nx*ny/(2^ngrids)):
global ngrids

% Define the initial guess 
guess = @(x,y) rand(length(x),length(y));

% Jacobi weighting:
w = 4/5;

%initialize a cell to hold the object on different grids
pA = cell(ngrids,1); 
pf = cell(ngrids,1);
pv = cell(ngrids,1);

%% Construct the initial guess and the RHS f:

% Define arrays of grid parameters
i = 0:ngrids-1;
nx = (nxfine+1)./(2.^i)-1;
ny = (nyfine+1)./(2.^i)-1;
dx = (bx-ax)./(nx+1);
dy = (by-ay)./(ny+1);

i = 1;
[pv{i}] = Build2DVector(guess,nx(i),ny(i),ax,bx,ay,by);
[pf{i}] = Build2DVector(f,nx(i),ny(i),ax,bx,ay,by);

%% Build the matrices on each level:
for i=1:ngrids
    [pA{i}] = BuildLaplacian2D(nx(i),ny(i),dx(i),dy(i));
end
%% Call Mu cycle:

% set up the exact solution so we can compare the error:
[Xmat,Ymat] = meshgrid(ax+dx(1):dx(1):bx-dx(1),ay+dy(1):dy(1):by-dy(1));
exact = pack_matToVec(exactV(Xmat,Ymat));
errornorm = zeros(nMucycles,1);
resnorm = zeros(nMucycles,1);
level = 1; %we start from level 1, the finest level

% Run the Mu cycles:
for i=1:nMucycles
    [pv] = MuCycle(pA,pv,pf,nx,nu1,nu2,w,level);
    % Compute the L2 norm after each MuCycle:
    errornorm(i) = sqrt(dx(1)^2*sum((exact-pv{level}).^2));
    resnorm(i) = sqrt(dx(1)^2*sum((pA{level}*pv{level}-pf{level}).^2));
end




