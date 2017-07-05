function [ lambda, z ] = InversePowerMethod ( A, kappa, epsilon, x, MaxIter )
%GORKA ERAÑA ROBLES - This function is an implementation of the inverse
%power iteration using the Rayleigh quotient method.
%   It follows the ideas developed in 2.2.

[m,~] = size(A);

Anorm = 1/norm(A,'fro');

for i = 1:MaxIter
    y = (A - kappa*eye(m))\x;
    ynorm = 1/norm(y);
    x1 = y*ynorm;
    w = x*ynorm;
    ro = (x1')*w;
    mu = kappa + ro;
    r = w - ro*x1;
    x = x1;
    kappa = mu;
    if (norm(r)*Anorm <= epsilon)
        lambda = mu;
        z = x;
        return
    end
end

error('Maximum number of iterations exceeded; increase option MaxIter.')

end

