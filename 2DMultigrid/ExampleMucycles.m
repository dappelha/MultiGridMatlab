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
% number of Mu cycles total to run:
nMucycles = 10;

% Problem size:
k = 5; % 9 gives 512x512 grids
% Specify the total number of interior fine points, (should be odd):
nx = 2^k -1 ; %best to just use nx=ny for now, there might be a bug.
ny = 2^k -1 ;
dx = (bx-ax)/(nx+1);
dy = dx;

% total number of grids:
global ngrids
ngrids = k;

% pv is a cell array that holds the solution (or error) at each grid level,
% pv{1} is the solution on the finest grid, i.e. the answer we are after
pv = cell(ngrids,1);

% initialize arrays that hold errors and error ratios (rho). 
errornorm = zeros(nMucycles,1);
resnorm = zeros(nMucycles,1);
resratio = zeros(nMucycles-1,1);
errorratio = zeros(nMucycles-1,1);
% Call the Multigrid routine which will run V or W cycles, not FMG:
[pv,resnorm,errornorm] = MultiGrid(ax,bx,ay,by,nx,ny,nMucycles,nu1,nu2,exactV,f);

%% Display the convergence:
display('residual norm = ')
fprintf(1,'%6.2E\n', resnorm)
resratio = resnorm(2:end)./resnorm(1:end-1);
display('r ratio = ')
fprintf(1,'%6.2f\n', resratio)
display('error norm = ')
fprintf(1,'%6.2E\n', errornorm)
errorratio = errornorm(2:end)./errornorm(1:end-1);
display('e ratio = ')
fprintf(1,'%6.2f\n', errorratio)

%% Plot the Solution and convergence:

i = 1;
xtabfull = ax:dx:bx;
ytabfull = ay:dy:by;
Vnewmat = zeros(nx+2,ny+2);
Vnewmat(2:end-1,2:end-1) = pack_vecToMat(pv{i},nx);
[Xmat,Ymat] = meshgrid(xtabfull,ytabfull);
figure
mesh(Xmat,Ymat,Vnewmat)
title(['After ', int2str(nMucycles), ' V(',int2str(nu1),',',int2str(nu2),') cycles'], 'Fontsize', 18);
set(gca,'FontSize', 14, 'FontName', 'Times New Roman');
zlabel('V')
xlabel('x')
ylabel('y')
% Plot the error and residual reduction:
figure
semilogy(1:nMucycles,errornorm)
title('Error Norm Reduction', 'Fontsize', 18);
set(gca,'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('V-cycle')
ylabel('||e||_2')
figure
semilogy(1:nMucycles,resnorm)
title('Residual Norm Reduction', 'Fontsize', 18);
set(gca,'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('V-cycle')
ylabel('||r||_2')
