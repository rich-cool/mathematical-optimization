%------------------Cholesky Factorization------------------%
function []=Cholesky(aa)

n = 5;
A = aa;
L = eye(n);
D = eye(n);

for j=1:n
    sum1 = 0;
    for s=1:j-1
        sum1 = sum1+(D(s,s)*L(j,s)^2);
    end
    D(j,j) = A(j,j) - sum1;
    for i = j+1:n
        sum2 = 0;
        for k=1:j-1
            sum2 = sum2+(D(k,k)*L(i,k)*L(j,k));
        end
        L(i,j) = (A(i,j)-sum2)/D(j,j);
    end
end

L
D
L'
A
L*D*L'
%-----------------------------------------------------------%