function [pv] = VisualMuCycle(pA,pv,pf,px,pAnalytic,pNumeric,nx,nu1,nu2,i)
%% Mu cycle scheme. i tells what level the scheme is solving on.
% Mu = 1 is a V cycle, Mu = 2 is a W cycle.
global ngrids 
global Mu
global vis


% If I'm not on the coarsest grid, relax and take the problem down a level
if i ~= ngrids
    if vis == 1
        pNumeric{i} = pA{i}\pf{i};
        figure
        plot(px{i},pAnalytic{i},'-',px{i},pNumeric{i},'.',px{i},pv{i},'x');
        title(['Exact error and approximate error before prerelaxation on grid ',num2str(i),'h'])
        legend('Exact Analytic Error','Numerical Error (direct solve)','Error Approximation (iterative)');
    end
    if vis == 2
        figure
        plot(px{i},pf{i}-pA{i}*pv{i},'--.');
        title(['Residual before prerelaxation on grid level ',num2str(i),'h'])
    end
    % Relax nu1 times on my current problem:
    % Jacobi weighting: 2/3 is ideal for 1-D.
    w = 2/3;
    [pv{i}] = WeightedJacobi(pv{i},pf{i},pA{i},nu1,w);
    if vis == 1
        figure
        plot(px{i},pAnalytic{i},'-',px{i},pNumeric{i},'.',px{i},pv{i},'x');
        title(['Exact error and approximate error after prerelaxation on grid ',num2str(i),'h'])
        legend('Exact Analytic Error','Numerical Error (direct solve)','Error Approximation (iterative)');
    end
    if vis == 2
        figure
        plot(px{i},pf{i}-pA{i}*pv{i},'--.');
        title(['Residual after prerelaxation on grid level ',num2str(i),'h'])
    end
    % compute residual and restrict it to coarse grid: (it replaces coarse f).
    pf{i+1} = Restrict1D( pf{i}-pA{i}*pv{i} ,nx(i));
    if vis == 2
        figure
        plot(px{i+1},pf{i+1},'--.');
        title([num2str(i),'h grid residual restricted to grid level ',num2str(i+1),'and used as RHS'])
    end
    % make the initual guess for the error on the coarser grid 0.
    pv{i+1} = zeros(nx(i+1),1);
    % compute the exact analytic error in our new approximation to the error
    % and restrict it to the coarser grid
    pAnalytic{i+1} = Restrict1D(pAnalytic{i}-pv{i},nx(i));
    % Call Mucycle scheme recursively mu times:
    for j=1:Mu
        pv = VisualMuCycle(pA,pv,pf,px,pAnalytic,pNumeric,nx,nu1,nu2,i+1);
    end
    % Prolongate back to the fine grid and add correction:
    pv{i} = pv{i} + Prolongate1D(pv{i+1},nx(i+1));
    if vis == 1
        figure
        plot(px{i},pAnalytic{i},px{i},pNumeric{i},'.',px{i},pv{i},'x');
        title(['Corrected solution (of error equation) on grid level ',num2str(i),'h, before postrelaxation'])
        legend('Exact Analytic Error','Numerical Error (direct solve)','Error Approximation (iterative)');
    end
    if vis == 2
        figure
        plot(px{i},pf{i}-pA{i}*pv{i},'--.');
        title(['Residual before post relaxation on grid level ',num2str(i),'h'])
    end
    % relax nu2 times after the correction:
    [pv{i}] = WeightedJacobi(pv{i},pf{i},pA{i},nu2,w);
    if vis == 1
        figure
        plot(px{i},pAnalytic{i},px{i},pNumeric{i},'.',px{i},pv{i},'x');
        title(['Corrected solution (of error equation) on grid level ',num2str(i),'h, after postrelaxation'])
        legend('Exact Analytic Error','Numerical Error (direct solve)','Error Approximation (iterative)');
    end
    if vis == 2
        figure
        plot(px{i},pf{i}-pA{i}*pv{i},'--.');
        title(['Residual after post relaxation on grid level ',num2str(i),'h'])
    end
else
    %if I'm already on the coarsest grid:

    % Solve thoroughly: (bunch of jacobi cycles are sufficient on very coarse grid)
        [pv{i}] = WeightedJacobi(pv{i},pf{i},pA{i},20,2/3);
%     pv{i} = pA{i}\pf{i}; %direct solve on coarsest grid is fast too,
%     and can be use instead of jacobi.
    if vis == 1
        pNumeric{i} = pA{i}\pf{i};
        figure
        plot(px{i},pAnalytic{i},px{i},pNumeric{i},'.',px{i},pv{i},'x');
        title(['Exact error and approximate error after solve on coarsest grid level (',num2str(i),'h)'])
        legend('Exact Analytic Error','Numerical Error (direct solve)','Error Approximation');
    end
end
