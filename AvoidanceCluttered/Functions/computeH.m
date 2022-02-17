function H = computeH(n, r)
    H = zeros(n-r+1);
    
    for i = 0:n-r
        for j = 0:n-r
            H(i+1, j+1) = (nchoosek(n-r, i)*nchoosek(n-r, j))/(nchoosek(2*n-r, i+j))/(2*n-r+1);
        end
    end
end