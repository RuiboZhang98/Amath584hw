function [Q,R] = MGramSchmidt(A)
% Modified Gram Schmidt
[m, n] = size(A); 
V = A;
Q = zeros(m,n); R = zeros(n,n);
    for i =1 : n
        v = V(:,i);
        r = sqrt(v'* v); %r_{ii} we use this to normalize v
        R(i,i) = r; 
        q = v./r;
        for j = i + 1 : n
            r = q' * V(:,j); %r_{ij} the projection of V(j) on q
            V(:, j) = V(:, j) - r * q; 
            R(i,j) = r;
        end
        Q(:, i) = q;
    end
end

