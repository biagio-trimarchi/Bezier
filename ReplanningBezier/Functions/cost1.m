function J = cost1(x, n, m, d, Pinit, Uinit)
    P = x(1:(n+1)*d);
    P = reshape(P, d, []);
    U = x((n+1)*d+1:end);
    
    J = 0;
    for k = 0:n
        Pdiff = P(:, k+1) - Pinit(:, k+1);
        J = J + Pdiff'*Pdiff;
    end
    for k = 0:m
        Udiff = U(k+1)-Uinit(k+1);
        J = J + Udiff^2;
    end
end