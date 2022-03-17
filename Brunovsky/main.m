%% %% %% %% BRUNOVSKY

clc
clear
addpath("./Functions")

n = 3;
m = 1;
s = 5;

D = computeD(s, n);
A = [0 1 0; 0 0 1; 0 0 0];
B = [0; 0; 1];
Y = [-1; -1; -1; 1; -10; 5];
U = s*(s-1)*(s-2)*D*Y;

x0 = [-1; 0; 0];

dt = 0.001;
t = 0:dt:1;

for tt = 1:length(t)
    y(tt) = 0;
    for k = 0:s
        y(tt) = y(tt) + bernsteinPol(s, k, (tt-1)*dt)*Y(k+1);
    end
    u(tt) = 0;
    for k = 0:s-n
        u(tt) = u(tt) + bernsteinPol(s-n, k, (tt-1)*dt)*U(k+1);
    end
end

x(:, 1) = x0;
for tt = 1:length(t)
    x(:, tt+1) = x(:, tt) + (A*x(:, tt) + B*u(tt))*dt;
end

plot(y)
hold on
%plot(u)
plot(x(1, :))
hold off

C = computeC()