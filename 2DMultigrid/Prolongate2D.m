function [vfine] = Prolongate2D(vcoarse,nx)
%% Takes a vector vcoarse and interpolates it to a fine grid vector vfine


% the n in the multigrid tutorial book is n = nx+1;
% the nx you send in here is the x dimension of vcoarse, do not confuse it
% with the larger x dimension of vfine.
% I use ghost zero points, so grow the vector 
ny = length(vcoarse)/nx;
vCoarseMat = zeros(nx+2,ny+2);
% 1st unpack the vector:
vCoarseMat(2:nx+1,2:ny+1) = pack_vecToMat(vcoarse,nx);
% Now allocate the fine grid 
vFineMat = zeros((nx+1)*2-1,(ny+1)*2-1);
% Do the interpolation:
vFineMat(2:2:end-1,2:2:end-1) = vCoarseMat(2:end-1,2:end-1);
vFineMat(1:2:end,2:2:end-1) = 0.5*(vCoarseMat(1:end-1,2:end-1) + vCoarseMat(2:end,2:end-1));
vFineMat(2:2:end,1:2:end) = 0.5*(vCoarseMat(2:end-1,1:end-1) + vCoarseMat(2:end-1,2:end));
vFineMat(1:2:end,1:2:end) = 0.25*(vCoarseMat(1:end-1,1:end-1)+ vCoarseMat(2:end,1:end-1)...
                                + vCoarseMat(1:end-1,2:end) + vCoarseMat(2:end,2:end));
% Pack the prolongated matrix back to a vector:
vfine = pack_matToVec(vFineMat);
