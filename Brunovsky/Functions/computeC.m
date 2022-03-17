function C = computeC(n, r)
    C = eye(n+1, n+1);
    
    for j = 1:r
        Caux = C;
        C = zeros(n+1-j, n+2-j);
        for k = 1:n+1-j
            C(k, :) = [ones(1, k), zeros(1, n+2-j - k)];
        end
        C = C*Caux;
    end
end