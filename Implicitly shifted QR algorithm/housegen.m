function [ u,nu ] = housegen( a )
%GORKA ERAÑA ROBLES - This function computes a vector u that generates a
%Householder reflection H = I - uu* satisfying Ha = nu·e1.
%   It follows the ideas developed in 3.4.1.
%   Accumulating this transformations any matrix can be reduced to upper
%   Hessenberg form.

u = a;
nu = norm(a);
if ( nu==0 )
    u(1) = sqrt(2);
    return
end
if ( u(1)~=0 )
    rho = (u(1)')/norm(u(1));
else
    rho = 1;
end

u = (rho/nu)*u;
u(1) = 1 + u(1);
u = u/sqrt(u(1));
nu = -(rho')*nu;

end