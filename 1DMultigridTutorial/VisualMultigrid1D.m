function [v] = VisualMultigrid1D(a,b,nxfine,nMucycles,nu1,nu2,analyticU,f)
%%

global vis
global ngrids

% Define the initial guess 
guess = @(x) rand(length(x),1);

%initialize a cell to hold the object on different grids
pA = cell(ngrids,1); %matrix operator
pf = cell(ngrids,1); % (RHS)
pv = cell(ngrids,1); % (unknowns)
px = cell(ngrids,1); %grid
pAnalytic = cell(ngrids,1); % the analytic solution
pNumeric = cell(ngrids,1); %The exact solution to the linear algebra problem

% Define arrays of gird parameters on each level
i = 0:ngrids-1;
nx = (nxfine+1)./(2.^i)-1;
dx = (b-a)./(nx+1);

% Build the matrices and grid on each level:
for i=1:ngrids
    % make the sparse second derivative matrix.
    ex = ones(nx(i),1); %array of ones
    pA{i} = spdiags([-ex 2*ex -ex], -1:1, nx(i), nx(i))/dx(i)^2;
    px{i} = a+dx(i):dx(i):b-dx(i); %only need interior grid points
end

% Initial approximation to the solution is just a starting guess (random).
v = guess(px{1}); 

% The analyticU is the analytic solution evaluated on the grid
% Similarly numericU is the gauss elim solution to the linear algebra problem:
numericU = pA{1}\f(px{1}); %SLOW, just for educational purpose.

if vis == 1 %Plot the initial guess and the exact numeric and analytic solutions
    figure
    plot(px{1},analyticU(px{1}),'-',px{1},numericU,'.',px{1},v,'x');
    title('Initial Iterative Approximation (random guess)')
    legend('Analytic Solution','Numerical Solution (direct solve)','Numerical Solution (initial guess)');
end


% Run the Mu cycles: (solving the error equation, Ae=r)

%To understand the process and the subsequent graphs, it is convenient to
%solve for the error instead of the solution directly. So pv holds the
%approximation for the error in the initial random guess.

residual = f(px{1})-pA{1}*v; %compute residual 
pv{1} = zeros(nx(1),1); %make 0 error guess (equivalent to random initial guess).
pf{1} = residual;
pAnalytic{1} = analyticU(px{1})-v; %the exact error difference with the analtyic solution

level = 1;
for i=1:nMucycles
    pv = VisualMuCycle(pA,pv,pf,px,pAnalytic,pNumeric,nx,nu1,nu2,level);
end
v = v + pv{1}; %correction to the initial random guess, u=v+e.
if vis == 1 %Plot the final solution and the exact solution
    figure
     plot(px{1},analyticU(px{1}),'-',px{1},numericU,'.',px{1},v,'x');
    title('Final Iterative Approximation After 1 Multigrid V-cycle')
    legend('Analytic Solution','Numerical Solution (direct solve)','Numerical Solution (Multigrid solve)');
end