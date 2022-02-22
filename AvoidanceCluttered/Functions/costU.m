function J = costU(x, m)
    J = 0;
    for i = 1:m
        J = J + (x(i+1)-x(i))^2;
    end
end