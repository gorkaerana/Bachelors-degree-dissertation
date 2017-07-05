function [ X ] = righteigvec( T )
%GORKA ERAÑA ROBLES - This routine computes, given an upper triangular 
%matrix T, its right eigenvectors. They are stored in the matrix X by 
%columns. They are normalized to have Frobenius norm one.
%   It follows the ideas developed in 4.5.
%   Input: upper triangular matrix T
%   Output: upper triangular matrix X of the right eigenvectors of T 

[n,~] = size(T);
smallnum = (n/eps)*realmin;
bignum = (eps/n)*realmax;
X = zeros(n);

for k =n:-1:1
    X(1:k-1,k) = -T(1:k-1,k);
    X(k,k) = 1;
    X(k+1:n,k) = 0;
    dmin = max(eps*abs(T(k,k)),smallnum);
    for j = k-1:-1:1
        d = T(j,j) - T(k,k);
        if ( abs(d) <= dmin )
            d = dmin;
        end
        if ( abs(X(j,k))/bignum >= abs(d) )
            s = abs(d)/abs(X(j,k));
            X(1:k,k) = s*X(1:k,k);
        end
        X(j,k) = X(j,k)/d;
        X(1:j-1,k) = X(1:j-1,k) - X(j,k)*T(1:j-1,j);
    end
    X(1:k,k) = X(1:k,k)/norm(X(1:k,k),'fro');
end

