% largely stolen from:
% https://uk.mathworks.com/matlabcentral/answers/368379-array-multiplication-along-certain-dimension
function B = matvec3dims(A,v)
%takes an A1xA2xN matrix A and multiplies the MxN matrix pages by an
%N-dimensional vector v
    [A1, A2, N] = size(A);
    D = permute(A,[3 1 2]);
    E = reshape(D,N,A1*A2);
    F = v.'*E;
    B = reshape(F,A1,A2).';
end

