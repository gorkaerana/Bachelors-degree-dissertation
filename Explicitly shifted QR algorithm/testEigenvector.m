% Script that computes the backwards error of the eigenvectors of a complex
% Schur form.

sim = 250;
c = 1;
for j = 1:sim
    n = randi([5,30]);
    A = exp(randn(n)*1i + randn(n));
    [T,Q] = complexschur(A);
    X = righteigvec(T);
    for i = 1:n
        semilogy(c,norm(T*X(:,i) - T(i,i)*X(:,i))/(norm(T)*norm(X(:,i))),'.');
        hold on
        c = c + 1;
    end
    
end