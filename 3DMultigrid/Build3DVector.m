function [ftab] = Build3DVector(f,nx,ny,nz,ax,bx,ay,by,az,bz)
%% Is used to construct a vector ftab from a 3D function handle f.

dx = (bx-ax)/(nx+1);
dy = (by-ay)/(ny+1);
dz = (bz-az)/(nz+1);
xtabfull = ax:dx:bx;
ytabfull = ay:dy:by;
ztabfull = az:dz:bz;
xtab = xtabfull(2:end-1);
ytab = ytabfull(2:end-1);
ztab = ztabfull(2:end-1);
% Define the initial guess and the RHS f:
[Xmat,Ymat,Zmat] = meshgrid(xtab,ytab,ztab);
fmat = f(Xmat,Ymat,Zmat);
ftab = reshape(fmat,[nx*ny*nz,1]);