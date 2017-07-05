function [ lambda, z ] = PowerMethod ( A, kappa, epsilon, x, MaxIter )
%GORKA ERAÑA ROBLES - This function is an implementation of the shifted
%power method.
%   It follows the ideas developed in 2.1.

Anorm = 1/norm(A,'fro');
xnorm = 1/norm(x);
x = x*xnorm;

for i = 1:MaxIter
    y = A*x;
    mu = (x')*y;
    r = y - mu*x;
    x = y - kappa*x;
    xnorm = 1/norm(x);
    x = x*xnorm;
    if (norm(r)*Anorm < epsilon)
        lambda = mu;
        z = x;
        return
    end
end

error('Maximum number of iterations exceeded; increase options MaxIter.');

end

