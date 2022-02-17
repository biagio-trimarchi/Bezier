function [c, ceq] = nonLinearConstraints2Agents(x, n, m, d, Q, Pd, Ud, Po, Uo, safeDist, A, b)
    % Useful quantities
    P1 = x(1:(n+1)*d);
    P1 = reshape(P1, d, []);
    U1 = x((n+1)*d+1:(n+1)*d+m+1);

    P2 = x((n+1)*d+m+2:(n+1)*d+m +(n+1)*d+1);
    P2 = reshape(P2, d, []);
    U2 = x((n+1)*d+m +(n+1)*d+2:end);

    % Equality Constraints
    
    ceq(1:d) = (P1(:, 2) - P1(:, 1))*(U1(2) - U1(1)) - (Pd(:, 2) - Pd(:, 1))*(Ud(2) - Ud(1));
    ceq(d+1:d+d) = (P1(:, end) - P1(:, end-1))*(U1(end) - U1(end-1)) - (Pd(:, end) - Pd(:, end-1))*(Ud(end) - Ud(end-1));
    ceq(d+d+1:d+d+d) = (P2(:, 2) - P2(:, 1))*(U2(2) - U2(1)) - (Po(:, 2) - Po(:, 1))*(Uo(2) - Uo(1));
    ceq(d+d+d+1:d+d+d+d) = (P2(:, end) - P2(:, end-1))*(U2(end) - U2(end-1)) - (Po(:, end) - Po(:, end-1))*(Uo(end) - Uo(end-1));
    %ceq = 0;
    
    % Inequality Constraints
    
    % Compute Composition P/U
    Cd1 = BezierComposition(P1(:, 1), P1(:, 2), P1(:, 3), P1(:, 4), P1(:, 5), P1(:, 6), U1(1), U1(2), U1(3), U1(4));
    Cd2 = BezierComposition(P2(:, 1), P2(:, 2), P2(:, 3), P2(:, 4), P2(:, 5), P2(:, 6), U2(1), U2(2), U2(3), U2(4));
    C1 = reshape(Cd1, d, []);
    C2 = reshape(Cd2, d, []);
    Cdiff = C1 - C2;
    
    % Establish Constraints
    c = zeros(2*n*m+1, 0);
    for k = 0:2*n*m
        c(k+1) = safeDist^2;
        for j = 1:d
            c(k+1) = c(k+1) - Cdiff(j, :)*Q(:, :, k+1)*Cdiff(j, :)';
        end
    end
    
    cDyn1 = A*Cd1 - b;
    cDyn2 = A*Cd2 - b;
%     % Velocity Constraint
%     auxV = 0;
%     for k = 0:n-1
%         auxV = auxV + norm(P(:, k+2)-P(:, k+1))^2;
%     end
%     auxU = 0;
%     for k = 0:m-1
%         auxU = auxU + norm(U(k+2) - U(k+1))^2;
%     end
%     
%     cV = auxV*auxU - vMax^2;
%     c(end) = cV;

    c = [c, cDyn1', cDyn2'];
    
end