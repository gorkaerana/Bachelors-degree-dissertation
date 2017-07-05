% Script that computes the backwards error of calculating the upper
% Hessenbger form of a matrix.

sim = 1000;
for j = 1:sim
    n = randi([5,30]);
    A = exp(randn(n)*1i + randn(n));
    [H,Q] = hessreduce(A);
    semilogy(j,norm(A - Q*H*Q')/norm(A),'.');
    hold on
end
axis([1 sim eps 10^(-14)])
ylabel('\epsilon_M < y < 10^{-14}')