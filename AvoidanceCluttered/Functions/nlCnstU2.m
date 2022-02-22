function [c, ceq] = nlCnstU2(x, P1, P2, n, m , d, dSafe, Q)
    ceq = 0;
    
    C1 = BezierComposition(P1(:, 1), P1(:, 2), P1(:, 3), P1(:, 4), P1(:, 5), P1(:, 6), x(1), x(2), x(3), x(4));
    C1 = reshape(C1, d, []);
    
    C2 = BezierComposition(P2(:, 1), P2(:, 2), P2(:, 3), P2(:, 4), P2(:, 5), P2(:, 6), x(5), x(6), x(7), x(8));
    C2 = reshape(C2, d, []);
    diff = C1 - C2;
    c = zeros(2*n*m+1, 1);
    for k = 0:2*n*m
        c(k+1) = dSafe^2;
        for j = 1:d
            c(k+1) = c(k+1) - diff(j, :)*Q(:, :, k+1)*diff(j, :)';
        end
    end
end