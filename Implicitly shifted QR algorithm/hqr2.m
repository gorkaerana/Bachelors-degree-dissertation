function [ H,Q ] = hqr2( H,Q,maxiter )
%GORKA ERAÑA ROBLES - This function applies the implicitly shifted QR 
%iteration to an upper Hessenberg form matrix H and applies the 
%transformations to the matrix Q.
%   It follows the ideas developed in 5.5.

[n,~] = size(H);
i2 = n;
iter = 0;

while i2 > 1
    
    if ( iter > maxiter ) % Return an error if the routine does not converge.
        error('Maximum number of iterations exceeded; increase option maxiter.')
    end
     
    oldi2 = i2;
    [H,Q,i1,i2] = backsearch2(H,Q,i2); % Find deflation rows.
    

    if ( i2 == oldi2 ) % If it does not converge, sum one to the counter.
        iter = iter + 1;
    else % If it does converge set the counter to 0.
        iter = 0;
    end
    
    % Apply the bulge chasing if the matrix hast not been completely
    % deflated.
    if ( i2 > 1 )
        u = startqr2step(H(i1,i1),H(i1,i1+1),H(i1+1,i1),H(i1+1,i1+1),H(i1+2,i1+1),H(i2-1,i2-1),H(i2-1,i2),H(i2,i2-1),H(i2,i2));
        [H,Q]= qr2step(H,Q,u,i1,i2);
    end

    
end

end

