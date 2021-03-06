function [c, ceq] = nlCnstP(x, Po1, Po2, n, d, dSafe, Q)
    ceq = 0;
    
    P1(1, :) = x(1:(n+1)); 
    P1(2, :) = x((n+1)*2+1:(n+1)*2+(n+1));
    
    diff1 = P1 - Po1;
    c1 = zeros(2*n+1, 1);
    for k = 0:2*n
        c1(k+1) = dSafe^2;
        for j = 1:d
            c1(k+1) = c1(k+1) - diff1(j, :)*Q(:, :, k+1)*diff1(j, :)';
        end
    end
    
    P2(1, :) = x((n+1)+1:2*(n+1));
    P2(2, :) = x((n+1)*2 +(n+1)+1:(n+1)*2 + 2*(n+1));
    
    diff2 = P2 - Po2;
    c2 = zeros(2*n+1, 1);
    for k = 0:2*n
        c2(k+1) = dSafe^2;
        for j = 1:d
            c2(k+1) = c2(k+1) - diff2(j, :)*Q(:, :, k+1)*diff2(j, :)';
        end
    end
    
    c = [c1; c2];
end