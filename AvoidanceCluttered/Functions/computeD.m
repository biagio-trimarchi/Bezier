function D = computeD(n, r)
    D = eye(n+1, n+1);
    
    for j = 1:r
        Daux = D;
        D = zeros(n+1-j, n+2-j);
        for k = 1:n+1-j
            D(k, :) = [zeros(1, k-1), -1, 1, zeros(1, n+2-j-2-(k-1))];
        end
        D = D*Daux;
    end
end