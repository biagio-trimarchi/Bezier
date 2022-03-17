function b = bernsteinPol(n, k, t)
    b = nchoosek(n, k)*t^k*(1-t)^(n-k);
end