function x = bezier(N, t, P, d)
    x = zeros(d, 1);
    for k = 0:N
        x = x + bernsteinPol(N, k, t)*P(:, k+1);
    end
end