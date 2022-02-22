function [c, ceq] = nlCnstUDyn(x, P, Co, n, m , d, dSafe, Q, Adyn, bdyn)
    ceq = 0;
    
    Cx = BezierComposition(P(:, 1), P(:, 2), P(:, 3), P(:, 4), P(:, 5), P(:, 6), x(1), x(2), x(3), x(4));
    C = reshape(Cx, d, []);
    diff = C - Co;
    c = zeros(2*n*m+1, 1);
    for k = 0:2*n*m
        c(k+1) = dSafe^2;
        for j = 1:d
            c(k+1) = c(k+1) - diff(j, :)*Q(:, :, k+1)*diff(j, :)';
        end
    end
    
    cdyn = Adyn*Cx - bdyn;
    
    c = [c; cdyn];
end