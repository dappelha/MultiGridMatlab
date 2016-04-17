%% Matlab version of Multigrid script 
%% solves the finite difference discretization of square 2-D Poisson.
clear all
close all
tic
% Visualize Data? vis = 1 for yes, vis = 0 for no.
vis = 1;

% Specify Simulation Parameters:

% Specify the domain:
ax = 0; bx = 1;
ay = 0; by = 1;
az = 0; bz = 1;
% Specify the total number of interior fine points, (should be odd):
nxfine = 63; %best to just use nx=ny for now, there might be a bug.
nyfine = nxfine;
nzfine = nxfine;
% Define the initial guess and the RHS functions:
kx = 5; ky = 5; kz = 5;
% guess = @(x,y,z)sin(kx*x*pi).*sin(ky*y*pi).*sin(kz*z*pi);
guess = @(x,y,z) rand(length(x),length(y),length(z));
% guess = @(x,y) sin(kx*x*pi).*sin(ky*y*pi) %+ sin(50*kx*x*pi).*sin(50*ky*y*pi);
% f = @(x,y) 2*((1-6*x.^2).*y.^2.*(1-y.^2)+(1-6*y.^2).*x.^2.*(1-x.^2));
% exactV = @(x,y) (x.^2 -x.^4).*(y.^4 - y.^2);
f = @(x,y,z) (sin(kx*x*pi).*sin(ky*y*pi).*sin(kz*z*pi));
exactV = @(x,y,z) (f(x,y,z)/(pi^2*kx^2 + pi^2*ky^2 + pi^2*kz^2));
%  f = @(x,y) 0*x+0*y;

% Specify the number of grids (must be <nx*ny*nz/(2^ngrids)):
global ngrids 
ngrids = 6;
% Jacobi weighting:
w = 6/7;
% Specify Mu, Mu=1 vcycle, Mu = 2 wcycle
global Mu
Mu = 1;
nMucycles = 10;
% Specify nu1 and nu2 for the number of prerelaxation and post relaxation
nu1 = 2; 
nu2 = 1;


%initialize a cell to hold the object on different grids
pA = cell(ngrids,1); 
pf = cell(ngrids,1);
pv = cell(ngrids,1);
presidual = cell(ngrids,1);

%% Construct the initial guess and the RHS f:

% Define arrays of grid parameters
i = 0:ngrids-1;
nx = (nxfine+1)./(2.^i)-1;
ny = (nyfine+1)./(2.^i)-1;
nz = (nzfine+1)./(2.^i)-1;
dx = (bx-ax)./(nx+1);
dy = (by-ay)./(ny+1);
dz = (bz-az)./(nz+1);
i = 1;
[pv{i}] = Build3DVector(guess,nx(i),ny(i),nz(i),ax,bx,ay,by,az,bz);
[pf{i}] = Build3DVector(f,nx(i),ny(i),nz(i),ax,bx,ay,by,az,bz);

% if vis == 1
%     %visualize the affects of restriction on a couple grids:
%     for i=1:ngrids % just a couple grids
%         % Assign the pointers:
%         [pv{i}] = Build3DVector(guess,nx(i),ny(i),nz(i),ax,bx,ay,by,az,bz);
%         [pf{i}] = Build3DVector(f,nx(i),ny(i),nz(i),ax,bx,ay,by,az,bz);
%         xtab = 0+dx(i):dx(i):1-dx(i);
%         ytab = 0+dy(i):dy(i):1-dy(i);
%         ztab = 0+dz(i):dz(i):1-dz(i);
%         [Xmat,Ymat,Zmat] = meshgrid(xtab,ytab,ztab);
%         figure
%         dummy = reshape(pv{i},[nx(i),ny(i),nz(i)]);
%         p = patch(isosurface(xtab,ytab,ztab,dummy,max(pv{i})/2));
% %         isonormals(xtab,ytab,ztab,dummy,p)
%         set(p,'FaceColor','red','EdgeColor','none');
%         daspect([1 1 1])
%         view(3); axis tight
%         camlight
% %         lighting gouraud
%     end
% end

%% Build the matrices on each level:
order = 2; % finite difference order
for i=1:ngrids
    [pA{i}] = fd3d(nx(i),ny(i),nz(i),order)/(dx(i)^2);
end
if vis==1
    for i=1:ngrids
        % Visualize the matrix:
        figure
        spy(pA{i})
    end
end
%% Call Mu cycle:
i = 1;
xtab = 0+dx(i):dx(i):1-dx(i);
ytab = 0+dy(i):dy(i):1-dy(i);
ztab = 0+dz(i):dz(i):1-dz(i);
[Xmat,Ymat,Zmat] = meshgrid(xtab,ytab,ztab);
exact = reshape(exactV(Xmat,Ymat,Zmat),[nx(i)*ny(i)*nz(i),1]);
errornorm = zeros(nMucycles,1);
resnorm = zeros(nMucycles,1);
level = 1; %we start from level 1
for i=1:nMucycles
    [pv] = MuCycle(pA,pv,pf,nx,ny,nz,nu1,nu2,w,level);
    % Compute the L2 norm after each MuCycle:
    errornorm(i) = sqrt(dx(1)^2*sum((exact-pv{level}).^2));
    resnorm(i) = sqrt(dx(1)^2*sum((pA{level}*pv{level}-pf{level}).^2));
end

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
toc
%% Plot the better guess:
if vis ==1
    i = 1;
    xtabfull = ax:dx(i):bx;
    ytabfull = ay:dy(i):by;
    ztabfull = az:dz(i):bz;
    Vnewmat = zeros(nxfine+2,nyfine+2,nzfine+2);
    Vnewmat(2:end-1,2:end-1,2:end-1) = reshape(pv{i},[nx(i),ny(i),nz(i)]);
    [Xmat,Ymat,Zmat] = meshgrid(xtabfull,ytabfull,ztabfull);
    
    figure
    IC = guess(Xmat,Ymat,Zmat);
    p = patch(isosurface(xtabfull,ytabfull,ztabfull,IC,max(max(max(IC)))/2));
    isonormals(xtabfull,ytabfull,ztabfull,IC,p)
    set(p,'FaceColor','red','EdgeColor','none');
    daspect([1 1 1])
    view(3); axis tight
    camlight
    lighting gouraud
    title('Initial Guess', 'Fontsize', 18);
    set(gca,'FontSize', 14, 'FontName', 'Times New Roman');
    zlabel('z')
    xlabel('x')
    ylabel('y')
    
    figure
    p = patch(isosurface(xtabfull,ytabfull,ztabfull,Vnewmat,max(max(max(Vnewmat)))/2));
    isonormals(xtabfull,ytabfull,ztabfull,Vnewmat,p)
    set(p,'FaceColor','red','EdgeColor','none');
    daspect([1 1 1])
    view(3); axis tight
    camlight
    lighting gouraud
    title(['After ', int2str(nMucycles), ' V(', int2str(nu1),',',int2str(nu2),') cycles'], 'Fontsize', 18);
    set(gca,'FontSize', 14, 'FontName', 'Times New Roman');
    zlabel('z')
    xlabel('x')
    ylabel('y')
    
    % Plot the error:
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
end
%% Direct solve with a zero RHS and 64x64 grid is 975 s
