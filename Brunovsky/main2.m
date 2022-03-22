%% %% %% %% BRUNOVSKY

clc
clear
addpath("./Functions")

n = 3;
m = 1;
s = 5;

a = [-1, -1, 3, 1];
A = [0 1 0; 0 0 1; -a(1:3)];
B = [0; 0; 1];
Y = [0; 1; 3; 1; -10; 5];

for i = 0:n
    E(i+1).E = computeE(s-i, i);
end

for i = 0:n
    D(i+1).D = computeD(s, i);
end

ED = zeros(s+1, s+1);
ED = ED + a(1)*E(1).E * D(1).D;
ED = ED + a(2)*E(2).E * D(2).D;
ED = ED + a(3)*E(3).E * D(3).D;
ED = ED + a(4)*E(4).E * D(4).D;

U = ED*Y;
Y1 = inv(ED)*U;

x0 = [0; 5; 20];

dt = 0.001;
t = 0:dt:1;

for tt = 1:length(t)
    y(tt) = 0;
    for k = 0:s
        y(tt) = y(tt) + bernsteinPol(s, k, (tt-1)*dt)*Y(k+1);
    end
    u(tt) = 0;
    for k = 0:s
        u(tt) = u(tt) + bernsteinPol(s, k, (tt-1)*dt)*U(k+1);
    end
end

x(:, 1) = x0;
for tt = 1:length(t)
    x(:, tt+1) = x(:, tt) + (A*x(:, tt) + B*u(tt))*dt;
end

plot(y)
hold on
%plot(u)
plot(real(x(1, :)))
hold off
legend('y', 'x(1)')