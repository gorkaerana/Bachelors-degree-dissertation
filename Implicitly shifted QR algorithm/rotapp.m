function [ x,y ] = rotapp( c,s,x,y )
%GORKA ERAÑA ROBLES - This function takes a plane rotation defined by c and
%s (the scalars returned by rotgen) and applies it to the vectors x and y.
%P·(x^t). It is implemented in complex arithmetic.
%  (y^t)
%   It follows the ideas developed in 4.4.
%   Input: rotation matrix P = (c   s ); vectors x and y
%                              (-s' c')
%   Output: x and y overwritten with P·(x^t)
%                                      (y^t)

t = c*x + s*y;
y = c*y - (s')*x;
x = t;

end

