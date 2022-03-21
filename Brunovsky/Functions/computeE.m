function E = computeE(n, r)
    E = zeros(n+r+1, n+1);
    
    for k = 0:n+r
        for j = max(0, k-r):min(n,k)
            E(k+1, j+1) = (nchoosek(r, k-j)* nchoosek(n, j))/nchoosek(n+r, k);
        end
    end
end