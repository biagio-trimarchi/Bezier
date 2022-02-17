%% Compute Bezier Composition
close all;
clear all;
clc;

% Curve Orders
m = 5;
n = 3;

P = sym('p', [m+1, 1]);
U = sym('u', [n+1, 1]);

A = P;

for kk = 1:m
    for ii = 0:m-kk
        for jj = 0:kk*n
            iIdx = ii+1;
            kIdx = kk+1;
            jIdx = jj+1;

            sum = 0;
            for ll = max([0, jj-n]):min([jj, kk*n-n])
                lIdx = ll+1;

                w1 = bin(kk*n-n, ll);
                w2 = bin(n, jj-ll);
                w3 = ((1-U(jIdx-lIdx+1,1))*A(iIdx, kIdx-1, lIdx)) + ...
                     (U(jIdx-lIdx+1,1)*A(iIdx+1, kIdx-1, lIdx));

                sum = sum + (w1*w2*w3);
            end

            A(iIdx, kIdx, jIdx) = (1/bin(kk*n, jj))*sum;
        end
    end
end

C(:,1) = A(1, m+1, 1:((m*n)+1));
C = simplify(C); % New Control Points

%% Test The Results -- Only For m = 5, n = 3
if m ~= 5 || n ~= 3
    return;
end

tic
% x Direction
p1 = 0;
p2 = 0;
p3 = 3;
p4 = 4;
p5 = 6;
p6 = 6;

u1 = 0;
u2 = 0.3;
u3 = 0.8;
u4 = 1.0;

cc_x = subs(C); % Evaluate Composed Control Points
pp_x = [p1;p2;p3;p4;p5;p6];

% y Direction
p1 = 0;
p2 = 0;
p3 = -2;
p4 = -1;
p5 = 0;
p6 = -2;

cc_y = subs(C);
pp_y = [p1;p2;p3;p4;p5;p6];

toc
% Unique Time Control Points (Can Be Different)
uu = [u1; u2; u3; u4];

r1 = [];
r2 = [];
for tt = 0:0.01:1
    for kk = 0:1:n
        bu(kk+1) = bin(n, kk)*(1-tt)^(n-kk)*tt^kk;
    end
    ut = bu*uu;

    for kk = 0:1:m
        bp1(kk+1) = bin(m, kk)*(1-ut)^(m-kk)*ut^kk;
    end
    r1 = cat(1, r1, bp1);

    for kk = 0:1:m*n
        bp2(kk+1) = bin(m*n, kk)*(1-tt)^(m*n-kk)*tt^kk;
    end
    r2 = cat(1, r2, bp2);
end

R1_x = r1*pp_x;
R1_y = r1*pp_y;
R2_x = r2*cc_x;
R2_y = r2*cc_y;

figure();
hold on;
plot(R1_x(:,1), R1_y(:,1), 'LineWidth', 1.5, 'Color', '#D95319');
plot(R2_x(:,1), R2_y(:,1), ':', 'LineWidth', 1.5, 'Color', '#EDB120');
plot(pp_x, pp_y, 'o', 'LineWidth', 1.5, 'Color', '#0072BD');
plot(pp_x, pp_y, 'LineWidth', 0.5, 'Color', '#0072BD');
plot(cc_x, cc_y, 'x', 'LineWidth', 1.5, 'Color', '#77AC30');
plot(cc_x, cc_y, 'LineWidth', 0.5, 'Color', '#77AC30');
hold off;

%% External Functions
function b = bin(x, y)
    b = factorial(x)/(factorial(y)*factorial(x-y));
end