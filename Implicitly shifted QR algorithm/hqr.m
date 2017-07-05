function [ H,Q ] = hqr( H,Q,maxiter )
%GORKA ERAÑA ROBLES - This routine overwrites H with a unitary similar 
%triangular matrix whose diagonals are the eigenvalues of H (the real Schur 
%form). Transformations are stored in Q.
%   It follows the ideas developed in 4.4.
%   Input: upper Hessenberg matrix H; number of maximum iterations maxiter
%   Output: Schur form overwritten in H; similarity transformation Q

[n,~] = size(H);
i2 = n;
iter = 0;
c = zeros(1,n);
s = zeros(1,n);

while 1
    iter = iter + 1;
    
    if (iter > maxiter) % Throws an error if maxiter is exceeded
        error('Maximum number of iterations exceeded; increase option maxiter.')
    end
    
    oldi2 = i2;
    [i1,i2] = backsearch(H,i2); % Check subdiagonal for near ceros, deflating points
    
    if ( i2==1 ) % End the function if H is upper triangular
        return
    end
    if ( i2~=oldi2 ) % Set iteration number to zero if there is another deflating row
        iter = 0;
    end
    
    % Compute Wilkinson shift
    kappa = wilkshift( H(i2-1,i2-1),H(i2-1,i2),H(i2,i2-1),H(i2,i2) );
    
    H(i1,i1) = H(i1,i1) - kappa; % Apply shift to the element of the diagonal that is left out of the loop
    for j = i1:i2-1 % Loop reducing the matrix to triangular form
        [ H(j,j),H(j+1,j),c(j),s(j) ] = rotgen( H(j,j),H(j+1,j) ); % Apply rotation so that the subdiagonal is set to zero
        H(j+1,j+1) = H(j+1,j+1) - kappa; % Apply shift to diagonal
        [ H(j,j+1:n),H(j+1,j+1:n) ] = rotapp( c(j),s(j),H(j,j+1:n),H(j+1,j+1:n) ); % Modify the involved rows
    end
    
    for k = i1:i2-1 % Loop applying the back multiplication
        [ H(1:k+1,k),H(1:k+1,k+1) ] = rotapp( c(k),conj(s(k)),H(1:k+1,k),H(1:k+1,k+1) );
        [ Q(1:n,k),Q(1:n,k+1) ] = rotapp( c(k),conj(s(k)),Q(1:n,k),Q(1:n,k+1) ); % Accumulate transformations
        H(k,k) = H(k,k) + kappa;
    end
    H(i2,i2) = H(i2,i2) + kappa; % 

end

