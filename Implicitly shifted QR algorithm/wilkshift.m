function [ kappa ] = wilkshift( a,b,c,d )
%GORKA ERAÑA ROBLES - This function computes the Wilkinson shift of a
%submatrix of order 2.
%   It is based on the ideas of 4.3.
%   Input: matrix B = ( a b )
%                     ( c d )
%   Output: shift kappa, nearest eigenvalue to d

kappa = d;
s = abs(a) + abs(b) + abs(c) + abs(d);

if (s == 0)
    return
end

q = (b/s)*(c/s);

if (q ~= 0)
    p = 0.5*((a/s) - (d/s));
    r = sqrt(p*p + q);
    if ( (real(p)*real(r) + imag(p)*imag(r)) < 0 )
        r = -r;
    end
    kappa = kappa - s*(q/(p+r));
end

end