function alpha = computeDistance(P, Q, alpha, tol)
upper = sqrt( ...
    max( [ ...
    norm(P(:, 1)-Q(:, 1)), ...
    norm(P(:, 1)-Q(:, end)), ...
    norm(P(:, end)-Q(:, 1)), ...
    norm(P(:, end)-Q(:, end)) ]));
if upper <= alpha
    alpha = upper;
end

if openGJK(P, Q) >= alpha*(1-tol)
    alpha
else
    [~, P1, P2] = castel(P, 1/2);
    [~, Q1, Q2] = castel(Q, 1/2);
    
    alpha = min([alpha, computeDistance(P1, Q1, alpha, tol)]);
    alpha = min([alpha, computeDistance(P1, Q2, alpha, tol)]);
    alpha = min([alpha, computeDistance(P2, Q1, alpha, tol)]);
    alpha = min([alpha, computeDistance(P2, Q2, alpha, tol)]);
end
end
