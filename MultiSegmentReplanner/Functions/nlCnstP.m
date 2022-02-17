function [c, ceq] = nlCnstP(x, Po,n, d, dSafe, Q)
    ceq = 0;
    
    P = [x(1:n+1); x((n+1)*d+1:(n+1)*(d+1))];
    
    diff = P - Po;
    c = zeros(2*n+1, 0);
    for k = 0:2*n
        c(k+1) = dSafe^2;
        for j = 1:d
            c(k+1) = c(k+1) - diff(j, :)*Q(:, :, k+1)*diff(j, :)';
        end
    end
    
end