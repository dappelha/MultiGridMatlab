function [A] = BuildLaplacian2D(nx,ny,dx,dy)
%% Returns the sparse matrix for 2D finite diffence laplacian:
Ix = speye(nx);
Iy = speye(ny);
ex = ones(nx,1);
ey = ones(ny,1);
% make the sparse second derivative matrix.
Bhx = spdiags([-ex 2*ex -ex], -1:1, nx, nx)/dx^2;
Bhy = spdiags([-ey 2*ey -ey], -1:1, ny, ny)/dy^2;
% Create the laplacian:
A = kron(Iy,Bhx) + kron(Bhy,Ix);
