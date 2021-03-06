function [c, ceq] = nlConstDyn2(x, P1, P2, U1, U2, Co1, n, m , d, dSafe, Q, Adyn, bdyn, T)
    ceq = n*m*(P1(:, end)-P1(:, end-1))*(1-x(1)) - (P2(:, 2)-P2(:, 1))*(x(2))/T;
   
    
    Cx1 = BezierComposition(P1(:, 1), P1(:, 2), P1(:, 3), P1(:, 4), P1(:, 5), P1(:, 6), 0, U1(2), x(1), 1);
    Cx2 = BezierComposition(P2(:, 1), P2(:, 2), P2(:, 3), P2(:, 4), P2(:, 5), P2(:, 6), 0, x(2), U2(3), 1);
    
    C1 = reshape(Cx1, d, []);
    diff = C1 - Co1;
    c1 = zeros(2*n*m+1, 1);
    for k = 0:2*n*m
        c1(k+1) = dSafe^2;
        for j = 1:d
            c1(k+1) = c1(k+1) - diff(j, :)*Q(:, :, k+1)*diff(j, :)';
        end
    end
    
%     cdyn1 = Adyn*Cx1 - bdyn;
%     cdyn2 = Adyn*Cx2 - bdyn;
    
%    c = [c1; cdyn1; cdyn2];
    c = c1;
end