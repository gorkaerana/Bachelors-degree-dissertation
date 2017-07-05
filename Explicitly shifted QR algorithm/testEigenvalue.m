% Script that computes the backwards error calculating the eigenvalues
% using the explicitly shifted QR iteration.

sim = 1000;
for j = 1:sim
    n = randi([5,30]);
    A = exp(randn(n)*1i + randn(n));
    [T,Q] = complexschur(A);
    eigenval = abs(eig(A));
    semilogy(j,norm(eigenval - abs(diag(T)))/norm(eigenval),'.');
    hold on
end
axis([1 sim eps 10^(-14)])
ylabel('\epsilon_M < y < 10^{-14}')