function [c, ceq] = nlCnstPU(x, Co,n, m, d, dSafe, Q, P0, U0)
    ceq = 0;
    
    P1 = [x(1:n+1); x((n+1)*d+1:(n+1)*(d+1))];
    U1 = x(2*d*(n+1)+1:(2*d*(n+1)+m+1));
    
    P2 = [x(n+2: 2*(n+1)); x((n+1)*(d+1)+1: (n+1)*(d+2))];
    U2 = x(2*d*(n+1)+m+2: end);
    
    ceq1 = (P1(:, 2) - P1(:, 1))*(U1(2) - U1(1)) - (P0(1:2, 2) - P0(1:2, 1))*(U0(1, 2) - U0(1, 1));     % initial velocity
    ceq2 = (P2(:, end) - P2(:, end-1))*(U1(end-1) - U1(end)) - (P0(3:4, end-1) - P0(3:4, end-1))*(U0(2, end) - U0(2, end-1));     % final velocity
    ceq3 = (P1(:, end) - P1(:, end-1))*(U1(end) - U1(end-1)) - (P2(:, 2) - P2(:, 1))*(U2(2) - U2(1));
    
    Cx = BezierComposition(P1(:, 1),P1(:, 2),P1(:, 3),P1(:, 4),P1(:, 5),P1(: ,6),U1(1),U1(2),U1(3),U1(4));
    C = reshape(Cx, d, []);
    
    diff = C - Co;
    c = zeros(2*(n*m)+1, 0);
    for k = 0:2*(n*m)
        c(k+1) = dSafe^2;
        for j = 1:d
            c(k+1) = c(k+1) - diff(j, :)*Q(:, :, k+1)*diff(j, :)';
        end
    end
    
    ceq = [ceq1, ceq2, ceq3];
end