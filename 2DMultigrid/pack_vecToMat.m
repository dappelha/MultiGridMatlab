function [Umat] = pack_vecToMat(Uvec,dim2)
%% Pack vector back to a matrix
% This function takes a vector that represents a tensor space and puts it
% back into matrix form. The dimensions of the desired matrix should be
% known beforehand, and the x dimension (dim2) should be sent to the
% function.
%% Now Smart, so can handle 2 or 3 coordinate data
if ndims(Uvec) == 2
    dim1 = size(Uvec,1)/dim2;
    dim3 = size(Uvec,2);
    Umat = zeros(dim1,dim2,dim3);
%     showme = Uvec
    for m = 1:dim1
        Umat(m, :,:) = Uvec((m-1)*dim2+1:m*dim2,:);
    end
end
if ndims(Uvec) == 1
    dim1 = length(Uvec)/dim2;
    Umat = zeros(dim1,dim2);
    for m = 1:dim1
        Umat(m, :) = Uvec((m-1)*dim2+1:m*dim2);
    end
end
% This way of unfolding corresponds to the description of folding given in
% pack_MatToVec.m, reproduced below:

% This way of folding takes x vectors at fixed y and stacks them in a
% vector. Say the x vector has 20 entries, then the first 20 entries in the
% packed vector will be the 20 x entries corresponding to the same fixed y
% value (the first y value). I call this type of packing a yx space. Then
% if you want to take the derivative wrt x you would use kron(I,B).