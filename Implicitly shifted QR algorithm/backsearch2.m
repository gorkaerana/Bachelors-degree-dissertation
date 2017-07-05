function [ H,Q,i1,i2 ] = backsearch2( H,Q,z )
%GORKA ERAÑA ROBLES - This function finds for deflating rows on a real
%Schur form matrix.
%   It is based on the ideas developed in 5.2.

normH = norm(H,'fro');
i1 = z;
i2 = z;

while i1>1

    if ( abs(H(i1,i1-1)) > eps*normH && i2 ~= 2 )
        i1 = i1 - 1; % Reduce i1
    else
        if i2 ~= 2
            H(i1,i1-1) = 0; % Deflate
        end

        % Check if it is a 1x1 or 2x2 block
        if ( i1 == i2 - 1 || i2 == 2 )
            
            % If it is a 2x2 block process it
            [H,Q] = blockprocess(H,Q,i2);
            
            if i2 ~= 2 % If it is a complex block go to row i1-1
               i2 = i1 - 1;
               i1 = i1 - 1; 
               
            else % If i2==2 then we have reached the firs 2x2 block
                i1 = 1;
                i2 = 1;
            end
            
        % If not, it is not a 2x2 block, it is a 1x1. Go to row i1-1 or
        % break the loop.
        elseif i1 == i2
            i2 = i1 - 1;
            i1 = i1 - 1; 
            else % Break the loop and finish the function
                break
        end
    end
end

end

