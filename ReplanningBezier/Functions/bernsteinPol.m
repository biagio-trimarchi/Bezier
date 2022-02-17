function b = bernsteinPol(n, k, t)
    b = binomial(n, k)*t^k*(1-t)^(n-k);
end