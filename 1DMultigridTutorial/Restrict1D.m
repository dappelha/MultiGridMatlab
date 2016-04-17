function [vcoarse] = Restrict1D(vfine,nx)
%% Takes a vector vfine and restricts it to a coarse grid vector vcoarse

% the n in the multigrid tutorial book is n = nx+1;

% Now allocate the coarse grid 
vcoarse = zeros((length(vfine)+1)/2-1,1);

% Do the restriction:
for i = 1:(nx+1)/2-1
    vcoarse(i) =  1/4*(vfine(2*i-1) + 2*vfine(2*i)+ vfine(2*i+1));
end
