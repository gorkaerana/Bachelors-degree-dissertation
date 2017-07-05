function [ H,Q ] = blockprocess( H,Q,i2 )
%GORKA ERAÑA ROBLES - This function processes a 2x2 block in backsearch2
%   If it is a 2x2 block with complex eigenvalues it does nothing, but if
%   it has real eigenvalues then it carries out a step of the QR iteration
%   on the block.

[n,~] = size(H);

sigma = wilkshift(H(i2-1,i2-1),H(i2-1,i2),H(i2,i2-1),H(i2,i2));

% Check if the eigenvalues are complex while the matrix is real
if ( ~isreal(sigma) && isreal([H(i2-1,i2-1),H(i2-1,i2);H(i2,i2-1),H(i2,i2)]) )
    return

% Process the block by applying a step of the QR routine on it.
else
    H = H - sigma*eye(n);
    [~,~,c,s] = rotgen(H(i2-1,i2-1)/norm(H(i2-1,i2-1)),H(i2,i2-1)/norm(H(i2,i2-1)));
%     [H(i2-1,i2-1),H(i2,i2-1),c,s] = rotgen(H(i2-1,i2-1),H(i2,i2-1));
    [H(i2-1,:),H(i2,:)] = rotapp(c,-s,H(i2-1,:),H(i2,:));
    [H(1:i2,i2-1),H(1:i2,i2)] = rotapp(c,-conj(s),H(1:i2,i2-1),H(1:i2,i2));
    [Q(:,i2-1),Q(:,i2)] = rotapp(c,-conj(s),Q(:,i2-1),Q(:,i2));
    H = H + sigma*eye(n);
    
end

end

