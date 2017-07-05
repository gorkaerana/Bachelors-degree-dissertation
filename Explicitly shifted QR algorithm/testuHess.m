% Script that computes how far is the transformation matrix Q from the
% reduction to Hessenberg form from being unitary.

sim = 1000;
for j = 1:sim
    n = randi([5,30]);
    A = exp(randn(n)*1i + randn(n));
    [H,Q] = hessreduce(A);
    semilogy(j,norm(eye(n) - Q*Q'),'.');
    semilogy(j,norm(eye(n) - Q'*Q),'o');
    hold on
end
axis([1 sim eps 10^(-14)])
ylabel('\epsilon_M < y < 10^{-14}')