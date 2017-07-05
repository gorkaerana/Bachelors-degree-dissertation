% Script that shows how close are the transformation matrices Q to being unitary.

sim = 1000;

for j = 1:sim
    
    n = randi([5 30]);
    A = randn(n);
    [~,Q] = realschur(A);
    semilogy(j,norm(eye(n) - Q*Q')/norm(A),'.');
    semilogy(j,norm(eye(n) - Q'*Q)/norm(A),'o');
    hold on
    
end

hold off