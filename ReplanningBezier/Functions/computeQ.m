function Q = computeQ(n)
    Q = zeros(n+1, n+1, 2*n+1);
    
    for k = 0:2*n
        for i = max(0, k-n):min(n,k)
            Q(k-i+1, i+1, k+1) = binomial(n, i)*binomial(n, k-i)/binomial(2*n, k);
        end
    end
end