function [ u ] = startqr2step( h11, h12, h21, h22, h32, hn1n1, hn1n, hnn1, hnn )
%GORKA ERAÑA ROBLES - This function returns the elements of the first
%column of H^2 - 2*Re(kappa)*H + |kappa|*I.
%   It is based on the ideas developed in 5.3.2.
%   Once it has computed the first three elements, it creates a vector u
%   that generates a Householder reflection so that the bulge chasing can
%   be applied to the rest of the matrix.

elements = [h11, h12, h21, h22, h32, hn1n1, hn1n, hnn1, hnn];
s = 1/max(abs(elements));
elements = s*elements;
p = elements(9) - elements(1);
q = elements(6) - elements(1);
r = elements(4) - elements(1);
c = [(p*q - elements(8)*elements(7))/elements(3) + elements(2); (r-p-q); elements(5)];
[u,~] = housegen(c);

end