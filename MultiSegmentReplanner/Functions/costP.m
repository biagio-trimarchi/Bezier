function [J, grad] = costP(x, n ,d, P0)
    P = [x(n+1); x(3*(n+1))];
    J = (P-P0)'*(P-P0);
    grad = P' - 2*P0';
end