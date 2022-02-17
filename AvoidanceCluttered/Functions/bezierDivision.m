function [P1, P2, Pdiv] = bezierDivision(P, t, n)
    Pdiv = zeros(n+1, n+1);
    
    Pdiv(1, :) = P;
    for k = 1:n
        for j = 0:n-k
            Pdiv(k+1, k+j+1) = (1-t)*Pdiv(k, k+j) + t*Pdiv(k, k+j+1);
        end
    end
    
    P1 = diag(Pdiv);
    P2 = flip(Pdiv(:, end));
end