function [c, ceq] = nlConstDyn2(x, P1, P2, U1, Co1, n, m , d, dSafe, Q, Adyn, bdyn)
    ceq(1) = (P1(end)-P1(end-1))*(1-x(1)) - (P2(2)-P2(1))*(x(2));
   
    
    Cx1 = BezierComposition(P1(:, 1), P1(:, 2), P1(:, 3), P1(:, 4), P1(:, 5), P1(:, 6), 0, U1(2), x(1), 1);
    
    C1 = reshape(Cx1, d, []);
    diff = C1 - Co1;
    c1 = zeros(2*n*m+1, 1);
    for k = 0:2*n*m
        c1(k+1) = dSafe^2;
        for j = 1:d
            c1(k+1) = c1(k+1) - diff(j, :)*Q(:, :, k+1)*diff(j, :)';
        end
    end
    
    c = c1;
    
end