function [ i1,i2 ] = backsearch( H,z )
%GORKA ERAÑA ROBLES - This function finds deflating rows on a complex Schur
%form matrix.
%   It is based on the ideas developed in 4.2.
%   Input: Hessenberg matrix H of order n; index z (1 < z <= n)
%   Output: indices i1 and i2 (i1, i2 <= z) holding one of the following 
%   conditions:
%        1) 1 <= i1 < i2 <= z, deflate at rows i1 and i2.
%        2) 1 = i2 = i2, matrix is completely deflated

i1 = z;
i2 = z;
normH = norm(H,'fro');

while (i1 > 1)
    
    if (abs(H(i1,i1-1)) < eps*normH) 
        H(i1,i1-1) = 0;
        if (i1 == i2)
            i2 = i1 - 1;
            i1 = i1 - 1;
        else
            return
        end
    else
        i1 = i1 - 1;
    end
    
end

end




