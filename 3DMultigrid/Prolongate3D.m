function [vfine] = Prolongate3D(vcoarse,nx,ny,nz)
%% Takes a vector vcoarse and interpolates it to a fine grid vector vfine


% the n in the multigrid tutorial book is n = nx+1;
% the nx you send in here is the x dimension of vcoarse, do not confuse it
% with the larger x dimension of vfine.
% I use ghost zero points, so grow the vector 
ny = length(vcoarse)/(nx*nz);
nz = length(vcoarse)/(nx*ny);
vCoarseMat = zeros(nx+2,ny+2,nz+2);
% 1st unpack the vector:
vCoarseMat(2:nx+1,2:ny+1,2:nz+1) = reshape(vcoarse,[nx,ny,nz]);
% Now allocate the fine grid 
vFineMat = zeros((nx+1)*2-1,(ny+1)*2-1, (nz+1)*2-1);
%% Do the interpolation:
% Corner points which correspond directly.
vFineMat(2:2:end-1,2:2:end-1,2:2:end-1) = vCoarseMat(2:end-1,2:end-1,2:end-1);
%% y midpoint of yz facets
vFineMat(2:2:end-1,1:2:end,2:2:end-1) = 0.5*(vCoarseMat(2:end-1,1:end-1,2:end-1) + vCoarseMat(2:end-1,2:end,2:end-1));
%% x midpoint of xz facets
vFineMat(1:2:end,2:2:end-1,2:2:end-1) = 0.5*(vCoarseMat(1:end-1,2:end-1,2:end-1) + vCoarseMat(2:end,2:end-1,2:end-1));
%% z midpoint of xz facets
vFineMat(2:2:end-1,2:2:end-1,1:2:end) = 0.5*(vCoarseMat(2:end-1,2:end-1,1:end-1) + vCoarseMat(2:end-1,2:end-1,2:end));

%% face midpoint of yz faces (constant x)
vFineMat(2:2:end-1,1:2:end,1:2:end) = 0.25*(vCoarseMat(2:end-1,1:end-1,1:end-1)+ vCoarseMat(2:end-1,2:end,1:end-1)...
                                + vCoarseMat(2:end-1,1:end-1,2:end) + vCoarseMat(2:end-1,2:end,2:end));

%% face midpoint of xz faces (constant y)
vFineMat(1:2:end,2:2:end-1,1:2:end) = 0.25*(vCoarseMat(1:end-1,2:end-1,1:end-1)+ vCoarseMat(2:end,2:end-1,1:end-1)...
                                + vCoarseMat(1:end-1,2:end-1,2:end) + vCoarseMat(2:end,2:end-1,2:end));
                           
%% face midpoint of xy faces (constant z)
vFineMat(1:2:end,1:2:end,2:2:end-1) = 0.25*(vCoarseMat(1:end-1,1:end-1,2:end-1)+ vCoarseMat(2:end,1:end-1,2:end-1)...
                                + vCoarseMat(1:end-1,2:end,2:end-1) + vCoarseMat(2:end,2:end,2:end-1));
                            
%% Center of cube:
vFineMat(1:2:end,1:2:end,1:2:end) = 1/8*(vCoarseMat(1:end-1,1:end-1,1:end-1) + vCoarseMat(2:end,1:end-1,1:end-1) + ...
                                         vCoarseMat(1:end-1,2:end,1:end-1) +   vCoarseMat(1:end-1,1:end-1,2:end) + ...
                                         vCoarseMat(1:end-1,2:end,2:end)   +   vCoarseMat(2:end,1:end-1,2:end) + ...
                                         vCoarseMat(2:end,2:end,1:end-1)   +   vCoarseMat(2:end,2:end,2:end));
%% Pack the prolongated matrix back to a vector:
vfine = reshape(vFineMat,[((nx+1)*2-1)*((ny+1)*2-1)*((nz+1)*2-1),1]);
