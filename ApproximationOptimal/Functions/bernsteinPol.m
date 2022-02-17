function b = bernsteinPol(N, k, t)
    b = binomial(N, k)*t^k*(1-t)^(N-k);
end