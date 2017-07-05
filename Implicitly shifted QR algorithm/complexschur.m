function [ T,Q ] = complexschur( A )
%GORKA ERA�A ROBLES - This function packs the reduction to Hessenberg form
%and the explicitly shifted QR routine in a singular function.

[H,Q] = hessreduce(A);
[T,Q] = hqr(H,Q,10000);

end

