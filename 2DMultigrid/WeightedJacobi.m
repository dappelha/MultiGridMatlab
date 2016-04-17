function [v] = WeightedJacobi(v,f,A,n_iterations,w)
%% Performs weighted jacobi iterations to solve Av=f.
D = diag(diag(A));
offD = A-D;
for i = 1:n_iterations
    v = (1-w)*v + w*(f - offD*v)./diag(A);
end