function [pv] = MuCycle(pA,pv,pf,nx,ny,nz,nu1,nu2,w,i)
%% Mu cycle scheme. i tells what level the scheme is solving on.
% Mu = 1 is a V cycle, Mu = 2 is a W cycle.
global ngrids 
global Mu
% Relax nu1 times on my current problem:
[pv{i}] = WeightedJacobi(pv{i},pf{i},pA{i},nu1,w);
% If I'm not on the coarsest grid, take the residual to a finer grid:
if i ~= ngrids
    % compute residual and restrict it to coarse grid: (it replaces coarse f).
    pf{i+1} = Restrict3D( pf{i}-pA{i}*pv{i} ,nx(i),ny(i),nz(i));
    % make the initual guess for the error on the coarser grid 0.
    pv{i+1} = zeros(nx(i+1)*ny(i+1)*nz(i+1),1);
    % Call Mucycle scheme recursively mu times:
    for j=1:Mu
        pv = MuCycle(pA,pv,pf,nx,ny,nz,nu1,nu2,w,i+1);
    end
    % Prolongate back to the fine grid and add correction:
    pv{i} = pv{i} + Prolongate3D(pv{i+1},nx(i+1),ny(i+1),nz(i+1));
else
    %if I'm already on the coarsest grid:
    pv{i} = pA{i}\pf{i};
    % Solve thoroughly: (20 jacobi cycles)
%     [pv{i}] = WeightedJacobi(pv{i},pf{i},pA{i},20,w);
end
% relax nu2 times on fine grid:
[pv{i}] = WeightedJacobi(pv{i},pf{i},pA{i},nu2,w);