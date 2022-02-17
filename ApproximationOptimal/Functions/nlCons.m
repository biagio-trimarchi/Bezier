function [c, ceq] = nlCons(xut, N, d, uMax, uMin, E, pO, nO, dP)
    ceq = 0;
    x = reshape(xut(1:(N+1)*d), d, []);
    u = reshape(xut((N+1)*d+1:end-1), d, []);
    tF = xut(end);
    
    xd = -N*diff(x,1,2);
    
    cf = zeros(1, N+1);
    cuu = zeros(1, N+1);
    cul = zeros(1, N+1);
    cE = zeros(1, nO*(N+1));
    
    for k = 0:N
        cf(k+1:k+2) = norm(bezier(N-1, k/N, xd, d) - tF*bezier(N, k/N, u, d)) - dP;
         cuu(k+1) = norm(bezier(N, k/N, u, d)) - uMax - dP;
         cul(k+1) = -norm(bezier(N, k/N, u, d)) + uMin + dP;
        for j = 1:nO
            cE(k*nO+j) = E - norm(bezier(N, k/N, x, d) - pO(:, j)) + dP;
        end
    end
    
    c = [cf, cuu, cul, cE];
end