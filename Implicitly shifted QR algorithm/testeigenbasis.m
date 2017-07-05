% Script that computes the backwards error of the eigenvectors of a real
% Schur form matrix.

simulations = 500;
counter = 1;
bad = 0;

for k = 1 : simulations
    n = randi([5 10]);
    A = randn(n);
    [T,~] = realschur(A);
    
    Y = eigenbasis(T);

    % subdiagonal is a vector representing the elements of the subdiagonal: 1
    % if it is nonzero, 0 if it is zero
    subdiagonal = not(not(diag(T,-1)));

    for j = 1:(n-1)

        if ( j ~= 1 && j ~= (n-1) )
            if ( subdiagonal(j+1) == 1 ) % The subdiagonal is in between two complex conjugate eigenvalues
                continue %  subdiagonal(j-1) == 1 &&
            elseif ( subdiagonal(j) == 1 )
                L = Y^(-1)*T*Y;
                a = norm(T*Y(:,j:j+1)-Y(:,j:j+1)*L(j:j+1,j:j+1))/(norm(T)*norm(Y(:,j:j+1)));
            else % Simple eigenvalue
                L = Y^(-1)*T*Y;
                a = norm(T*Y(:,j+1)-Y(:,j+1)*L(j+1,j+1))/(norm(T)*norm(Y(:,j+1)));
            end

        elseif ( j == 1 )
            if ( subdiagonal(1) == 1 ) % Complex eigenvalue
                L = Y^(-1)*T*Y;
                a = norm(T*Y(:,1:2)-Y(:,1:2)*L(1:2,1:2))/(norm(T)*norm(Y(:,1:2)));
            else % Simple eigenvalue
                L = Y^(-1)*T*Y;
                a = norm(T*Y(:,1)-Y(:,1)*L(1,1))/(norm(T)*norm(Y(:,1)));
            end

        elseif ( j == (n-1) )
            if ( subdiagonal(n-1) == 1 ) % Complex eigenvalue
                L = Y^(-1)*T*Y;
                a = norm(T*Y(:,n-1:n)-Y(:,n-1:n)*L(n-1:n,n-1:n))/(norm(T)*norm(Y(:,n-1:n)));
            else % Simple eigenvalue
                L = Y^(-1)*T*Y;
                a = norm(T*Y(:,n)-Y(:,n)*L(n,n))/(norm(T)*norm(Y(:,n)));
            end

        end

        if (a > 10^(-13))
            bad = bad + 1;
            disp(A);
            return
        end

        semilogy(counter,a,'.')
        hold on
        counter = counter + 1;

    end
        
end

ylabel('10^{-18} < y < 10^{-15}')
% axis([1 counter 10^(-18) 10^(-15)])

hold off
