function [U,S,V] = RSVD(A,k)
%   randomized SVD
%   k is the sampling constant.
Omega = rand(size(A,2),k);
[Q,~] = qr(A*Omega,0);
B = Q'*A;
[OU,S,V] = svd(B,'econ');
U = Q * OU;
end

