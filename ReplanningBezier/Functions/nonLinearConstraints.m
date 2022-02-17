function [c, ceq] = nonLinearConstraints(x, n, m, d, Q, Co, Pd, Ud, safeDist, A, b)
    % Useful quantities
    P = x(1:(n+1)*d);
    P = reshape(P, d, []);
    U = x((n+1)*d+1:end);

    % Equality Constraints
    
    ceq(1:d) = (P(:, 2) - P(:, 1))*(U(2) - U(1)) - (Pd(:, 2) - Pd(:, 1))*(Ud(2) - Ud(1));
    %ceq(d+1:d+d) = (P(:, end) - P(:, end-1))*(U(end) - U(end-1)) - (Pd(:, end) - Pd(:, end-1))*(Ud(end) - Ud(end-1));
    %ceq1 = (P(:, 2) - P(:, 1))*(U(2) - U(1)) - (Cd(:, 2) - Cd(:, 1));
    %ceq2 = (P(:, end) - P(:, end-1))*(U(end) - U(end-1));
    
    %ceq = [ceq1, ceq2];
    % Inequality Constraints
    
    % Compute Composition P/U
    
    C = BezierComposition33(P(:, 1), P(:, 2), P(:, 3), P(:, 4), U(1), U(2), U(3), U(4));
    Cx = reshape(C, [], 1);
    Cdiff = C - Co;
    
    % Establish Constraints
    c = zeros(2*n*m+1, 0);
    for k = 0:2*n*m
        c(k+1) = safeDist^2;
        for j = 1:d
            c(k+1) = c(k+1) - Cdiff(j, :)*Q(:, :, k+1)*Cdiff(j, :)';
        end
    end
    
    cDyn = A*Cx - b;
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

    c = [c, cDyn'];
    
end