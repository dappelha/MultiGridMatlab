function [vcoarse] = Restrict2D(vfine,nx)
%% Takes a vector vfine and restricts it to a coarse grid vector vcoarse

% the n in the multigrid tutorial book is n = nx+1;
% 1st unpack the vector:
vFineMat = pack_vecToMat(vfine,nx);
% Now allocate the coarse grid 
ny = length(vfine)/nx;
vCoarseMat = zeros((size(vFineMat)+1)/2-1);

% Do the restriction:
for i = 1:(nx+1)/2-1
    for j=1:(ny+1)/2-1
        vCoarseMat(i,j) =  1/16*(vFineMat(2*i-1,2*j-1) + vFineMat(2*i-1,2*j+1)...
            + vFineMat(2*i+1,2*j-1) + vFineMat(2*i+1,2*j+1) ...
            + 2*(vFineMat(2*i,2*j-1) + vFineMat(2*i,2*j+1) ...
            +    vFineMat(2*i-1,2*j) + vFineMat(2*i+1,2*j))...
            + 4*vFineMat(2*i,2*j));
    end
end
% Pack the restricted matrix back to a vector:
vcoarse = pack_matToVec(vCoarseMat);