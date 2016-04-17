close all
clear all
% Specify Simulation Parameters:

% Specify the domain:
ax = 0; bx = 1;
ay = 0; by = 1;

% Define the RHS function and the known analytic solution:
f = @(x,y) 2*((1-6*x.^2).*y.^2.*(1-y.^2)+(1-6*y.^2).*x.^2.*(1-x.^2));
exactV = @(x,y) (x.^2 -x.^4).*(y.^4 - y.^2);

% Specify Mu, Mu=1 vcycle, Mu = 2 wcycle
global Mu
Mu = 1;
% Specify the number of pre and post relaxations i.e. V(nu1,nu2).
nu1 = 2;
nu2 = 2;
% number of Mu cycles at each grid level in FMG:
ncycles = 1;


global ngrids
k = 5; % 9 gives 512x512 grids

% Specify the total number of interior fine points, (should be odd):
nx = 2^k -1 ; %best to just use nx=ny for now, there might be a bug.
ny = 2^k -1 ;
dx = (bx-ax)/(nx+1);
dy = dx;
ngrids = k;
% pv is a cell array that holds the solution (or error) at each grid level,
% pv{1} is the solution on the finest grid, i.e. the answer we are after
pv = cell(ngrids,1);

% Call FullMultigrid routine:
pv = runFMG(ax,bx,ay,by,nx,ny,ncycles,nu1,nu2,f);
% Compute the L2 error:
[Xmat,Ymat] = meshgrid(ax+dx:dx:bx-dx,ay+dy:dy:by-dy);
exact = pack_matToVec(exactV(Xmat,Ymat));

L2error = sqrt(dx^2*sum((exact-pv{1}).^2));
display(['The L2 error is ', num2str(L2error)])


%% Plot the Solution and pointwise error:

i = 1;
xtabfull = ax:dx:bx;
ytabfull = ay:dy:by;
Vnewmat = zeros(nx+2,ny+2);
Vnewmat(2:end-1,2:end-1) = pack_vecToMat(pv{i},nx);
[Xmat,Ymat] = meshgrid(xtabfull,ytabfull);
figure
mesh(Xmat,Ymat,Vnewmat)
title(['Solution after FMG with', ' V(',int2str(nu1),',',int2str(nu2),') cycles'], 'Fontsize', 18);
set(gca,'FontSize', 14, 'FontName', 'Times New Roman');
zlabel('V')
xlabel('x')
ylabel('y')
% Plot the pointwise error:
figure
errorMat = zeros(nx+2,ny+2);
errorMat(2:end-1,2:end-1) = pack_vecToMat(abs(pv{i}-exact),nx(i));
mesh(Xmat,Ymat,errorMat)
title(['Pointwise error after FMG with', ' V(',int2str(nu1),',',int2str(nu2),') cycles'], 'Fontsize', 18);
zlabel('|u_{*}-u_{fmg}|')
xlabel('x')
ylabel('y')

