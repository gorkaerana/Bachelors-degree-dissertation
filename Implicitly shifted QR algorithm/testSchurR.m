% Script that computes the backwards error of compouting the real Schur
% form using the implicitly shifted QR algorithm.

sim = 1000;

for j = 1:sim
    
    n = randi([5 30]);
    A = randn(n);
    [T,Q] = realschur(A);
    nn = norm(A - Q*T*Q')/norm(A);
    semilogy(j,nn,'.');
    hold on
    
end

hold off