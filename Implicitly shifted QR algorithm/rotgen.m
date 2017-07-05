function [ a,b,c,s ] = rotgen( a,b )
%GORKA ERAÑA ROBLES - This function generates a Givens rotation from
%elements a and b. It is implemented in complex arithmetic.
%   It follows the ideas developed in 4.4.
%   Input: quantities a, b where (c   s )·(a) = (nu*a/abs(a))
%                                (-s' c') (b)   (     0     )
%   Output: constants c and s; overwrites a with its final version and b with
%   zero

if ( b==0 )
    c = 1;
    s = 0;
    return
end
if ( a==0 )
    c = 0;
    s = 1;
    a = b;
    b = 0;
    return
end

mu = a/abs(a);
tau = abs(real(a)) + abs(imag(a)) + abs(real(b)) + abs(imag(b));
nu = tau*sqrt(abs(a/tau)^2 + abs(b/tau)^2);
c = abs(a)/nu;
s = mu*(b')/nu;
a = nu*mu;
b = 0;

end

