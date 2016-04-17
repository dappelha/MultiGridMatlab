function [pv] = FullMultiGrid(pA,pv,pf,nx,ncycles,nu1,nu2,w,i)
%% Author David Appelhans 2011
global ngrids
% if I'm not on the coarsest grid
if i~= ngrids
    pv = FullMultiGrid(pA,pv,pf,nx,ncycles,nu1,nu2,w,i+1);
    pv{i} = Prolongate2D(pv{i+1},nx(i+1));
else
    % If i'm on the coarsest grid just guess whatever and then V-cycle
    pv{i} = rand(nx(i)^2,1);
end
% Solve on the current grid
for j = 1:ncycles
    pv = MuCycle(pA,pv,pf,nx,nu1,nu2,w,i);
end

    