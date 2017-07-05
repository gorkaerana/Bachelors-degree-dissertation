function [ T,Q ] = realschur( A )
%GORKA ERAÑA ROBLES - This function packs the Hessenberg reduction and the
%explicitly shifted QR iteration into a singular function.

[H,Q] = hessreduce(A);
[T,Q] = hqr2(H,Q,10000);

end

