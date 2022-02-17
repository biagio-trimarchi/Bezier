function J = cost2Agents(x, n, m, d, Pinit1, Uinit1, Pinit2, Uinit2)
P1 = x(1:(n+1)*d);
P1 = reshape(P1, d, []);
U1 = x((n+1)*d+1:(n+1)*d+m+1);

P2 = x((n+1)*d+m+2:(n+1)*d+m +(n+1)*d+1);
P2 = reshape(P2, d, []);
U2 = x((n+1)*d+m +(n+1)*d+2:end);

J = 0;
for k = 0:n
    Pdiff = P1(:, k+1) - Pinit1(:, k+1);
    J = J + Pdiff'*Pdiff;
end
% for k = 0:m
%     Udiff = U1(k+1)-Uinit1(k+1);
%     J = J + Udiff^2;
% end
for k = 0:n
    Pdiff = P2(:, k+1) - Pinit2(:, k+1);
    J = J + Pdiff'*Pdiff;
end
% for k = 0:m
%     Udiff = U2(k+1)-Uinit2(k+1);
%     J = J + Udiff^2;
% end
end