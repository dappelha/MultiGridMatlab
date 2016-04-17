function [ftab] = Build2DVector(f,nx,ny,ax,bx,ay,by)
%% Is used to construct a vector ftab from a function 2D function handle f.
dx = (bx-ax)/(nx+1);
dy = (by-ay)/(ny+1);
xtabfull = ax:dx:bx;
ytabfull = ay:dy:by;
xtab = xtabfull(2:end-1);
ytab = ytabfull(2:end-1);
% Define the initial guess and the RHS f:
[Xmat,Ymat] = meshgrid(xtab,ytab);
fmat = f(Xmat,Ymat);
ftab = pack_matToVec(fmat);