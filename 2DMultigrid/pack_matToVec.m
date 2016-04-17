function [Uvec] = pack_matToVec(Umat)
%% Pack matrix into a vector
% This function takes a matrix with rows and columns corresponding to 
% y and x axis respectively, and packs it into a vector. If the matrix 
% is N by N, then the vector will have N^2 entries.
%% Now Smart, so can handle 2 or 3 coordinate data
if ndims(Umat) == 3
    dims = size(Umat);
    Uvec = zeros(dims(1)*dims(2), dims(3));
    for m = 1:dims(1)
        Uvec((m-1)*dims(2)+1:m*dims(2),:) = Umat(m, :,:);
    end
end

if ndims(Umat) == 2
    dims = size(Umat);
    Uvec = zeros(dims(1)*dims(2), 1);
    for m = 1:dims(1)
        Uvec((m-1)*dims(2)+1:m*dims(2)) = Umat(m, :);
    end
end
% This way of folding takes x vectors at fixed y and stacks them in a
% vector. Say the x vector has 20 entries, then the first 20 entries in the
% packed vector will be the 20 x entries corresponding to the same fixed y
% value (the first y value). I call this type of packing a yx space. Then
% if you want to take the derivative wrt x you would use kron(I,B).