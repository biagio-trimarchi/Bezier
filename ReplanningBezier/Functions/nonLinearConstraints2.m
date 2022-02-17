function [c, ceq] = nonLinearConstraints2(x, n, d, Co, Q, safeDist)
    % Useful quantities
    C = reshape(x, d, []);
    % Inequality Constraints
    Cdiff = C - Co;
    
    % Establish Constraints
    c = zeros(2*n+1, 0);
    for k = 0:2*n
        c(k+1) = safeDist^2;
        for j = 1:d
            c(k+1) = c(k+1) - Cdiff(j, :)*Q(:, :, k+1)*Cdiff(j, :)';
        end
    end
    
    ceq = 0;
end