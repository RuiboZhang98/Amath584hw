function [lambda,v] = PowerIter(A,v0,iter)
% power iteration method of matrix A
% iter: the number of iterations. v0:initial guess.
% lambda: eigenvalue v:unit length eigenvector
v = v0;
v = v/norm(v);
for k = 1 : iter
    w = A*v;
    v = w/norm(w);
    lambda = v'*A*v;
end
end

