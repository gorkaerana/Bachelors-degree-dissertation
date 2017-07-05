function [ Y ] = eigenbasis( T )
%GORKA ERAÑA ROBLES - This function computes the eigenbases of a real Schur
%form matrix T.
%   The implementation follows the computations developed in section 5.6 of
%   the dissertation.

[n,~] = size(T);
j = n;
Y = zeros(n);

while ( j ~= 0 )
        
    if ( j == 1 ) % Base case number 1: the first eigenvalues is simple
        Y(1,1) = 1;
        j = j - 1;
    
    elseif ( T(2,1) ~= 0 && j == 2 ) % Base case number 2: the first two eigenvalues are complex conjugate
        [L,~] = complexschur(T(1:2,1:2));
        L = [real(L(1,1)),imag(L(1,1));-imag(L(1,1)),real(L(1,1))];
        Y(1:2,1:2) = L;
%         Y(1:2,1:2) = eye(2);
        j = j - 2;
    
    elseif ( T(j,j-1) == 0 && j > 1 ) % The current eigenvalue is simple
        
        if ( j == 2 )
            lambda = T(j,j);
            Y(j,j) = 1;
            Y(j-1,j) = T(j-1,j)/(lambda - T(j-1,j-1));
            Y(1:j,j) = Y(1:j,j)/norm(Y(1:j,j),'fro');
            j = j - 1;
            
        elseif ( T(j-1,j-2) == 0 && j > 2 ) % The eigenvalue above is simple
            lambda = T(j,j);
            Y(j,j) = 1;
            Y(j-1,j) = T(j-1,j)/(lambda - T(j-1,j-1));
            Y(1:j-2,j) = (T(1:j-2,1:j-2)-lambda*eye(j-2))\(-T(1:j-2,j-1)*Y(j-1,j)-T(1:j-2,j));
            Y(1:j,j) = Y(1:j,j)/norm(Y(1:j,j),'fro');
            j = j - 1;
            
        elseif ( j >= 3 ) % The eigenvalues above are complex conjugate
            lambda = T(j,j);
            Y(j,j) = 1;
            Y(j-2:j-1,j) = (T(j-2:j-1,j-2:j-1)-lambda*eye(2))\(-T(j-2:j-1,j));
            if ( j > 3 )
                Y(1:j-3,j) = (T(1:j-3,1:j-3)-lambda*eye(j-3))\(-T(1:j-3,j-2:j-1)*Y(j-2:j-1,j)-T(1:j-3,j));
            end
            Y(1:j,j) = Y(1:j,j)/norm(Y(1:j,j),'fro');
            j = j - 1;
            
        end
        
    else % The current eigenvalues are complex conjugates
        
        if ( j == 3 ) % Base case: the one above is the last simple eigenvalue
            Y(1,1) = 1;
            [L,Z] = complexschur(T(j-1:j,j-1:j));
            X3 = righteigvec(L);
            X3 = Z*X3;
            L = [real(L(1,1)),imag(L(1,1));-imag(L(1,1)),real(L(1,1))];
            Y(2:j,2:j) = [real(X3(:,1)),imag(X3(:,1))];
            Y(1,2:j) = T(1,2:j)*Y(2:j,2:j)/(L-T(1,1)*eye(2));
            Y(:,j-1:j) = Y(:,j-1:j)/norm(Y(:,j-1:j),'fro');
            j = j - 3;
            
        elseif ( T(j-2,j-3) == 0 && j > 3 ) % The eigenvalue above is simple
            [L,Q1] = complexschur(T(j-1:j,j-1:j));
            X3 = righteigvec(L);
            X3 = Q1*X3;
            X3 = [real(X3(:,1)),imag(X3(:,1))]; % The definitive X3
            Y(j-1:j,j-1:j) = X3;
            L = [real(L(1,1)),imag(L(1,1));-imag(L(1,1)),real(L(1,1))]; % Instead of computing X3^(-1)*T33*X3, we follow theoretical results
            Y(j-2,j-1:j) = (T(j-2,j-1:j)*X3)/(L - T(j-2,j-2)*eye(2));
            x2 = Y(j-2,j-1:j);
            Y(1:j-3,j-1:j) = sylvester(T(1:j-3,1:j-3),-L,-T(1:j-3,j-1:j)*X3-T(1:j-3,j-2)*x2);
            Y(:,j-1:j) = Y(:,j-1:j)/norm(Y(:,j-1:j),'fro');
            j = j - 2;
            
        elseif ( j >= 4 ) % The eigenvalues above are complex conjugates
            if ( j > 4 )
                [L,Q1] = complexschur(T(j-1:j,j-1:j));
                X3 = righteigvec(L);
                X3 = Q1*X3;
                X3 = [real(X3(:,1)),imag(X3(:,1))]; % The definitive X3
                Y(j-1:j,j-1:j) = X3;
                L = [real(L(1,1)),imag(L(1,1));-imag(L(1,1)),real(L(1,1))]; % Instead of computing X3^(-1)*T33*X3, we follow theoretical results
                Y(j-3:j-2,j-1:j) = sylvester(T(j-3:j-2,j-3:j-2),-L,-T(j-3:j-2,j-1:j)*X3);
                X2 = Y(j-3:j-2,j-1:j);
                Y(1:j-4,j-1:j) = sylvester(T(1:j-4,1:j-4),-L,-T(1:j-4,j-3:j-2)*X2-T(1:j-4,j-1:j)*X3);
            else % Base case: both eigenvalues above are complex conjugates
                Y(3:4,3:4) = eye(2);
                Y(1:2,3:4) = sylvester(T(1:2,1:2),-T(3:4,3:4),-T(1:2,3:4));
            end
            Y(:,j-1:j) = Y(:,j-1:j)/norm(Y(:,j-1:j),'fro');
            j = j - 2;
            
        end
        
    end
    
end

end

