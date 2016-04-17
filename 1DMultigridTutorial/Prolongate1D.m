function [vfine] = Prolongate1D(vcoarse,nx)
%% Takes a vector vcoarse and interpolates it to a fine grid vector vfine


% the n in the multigrid tutorial book is n = nx+1;
% the nx you send in here is the x dimension of vcoarse, do not confuse it
% with the larger x dimension of vfine.
% I use ghost zero points, so grow the vector 
vCoarse = zeros(1,nx+2);
% 1st unpack the vector:
vCoarse(2:nx+1) = vcoarse;
% Now allocate the fine grid 
vfine = zeros((nx+1)*2-1,1);
% Do the interpolation:
vfine(2:2:end-1) = vCoarse(2:end-1);
vfine(1:2:end) = 0.5*(vCoarse(1:end-1) + vCoarse(2:end));