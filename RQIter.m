function [lambda,v] = RQIter(A,v0,iter,tol)
%   Rayleigh Quotient iteration for matrix A
%   v0:initial guess. eps:convergent constant. 
%   lambda:eigenvalue.v:eigenvector
m = size(A,1);
v = v0/norm(v0);
lambda = v'*A*v;
for i = 1:iter
    w = (A - lambda*eye(m))\v;
    v1 = v;
    v = w/norm(w);
    lambda = v'*A*v;
    err = norm(v-v1);
    i = i+1;
    if err<tol,break;end
end
end

