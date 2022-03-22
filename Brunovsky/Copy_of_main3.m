%% %% %% %% BRUNOVSKY

clc
clear
addpath("./Functions")

%% PARAMETERS
r = 3;
n = 5;
s = 7;

%% MATRICES
a4 = [1, 1, 10, 1];
a3 = [-30, 10];
A1 = [0, 1; -1, -1];
A2 = [0, 0, 0; 0, 0, 0];
A3 = [0, 0; 0, 0; a3];
A4 = [0 1 0; 0 0 1; -a4(1:3)];
B = [0; 0; 1];
Y = [0; 0; 0; 3; 1; -10; -10; 5];

x0 = [0; 0; 0];
z0 = [1; 10];

A = [A1, A2; A3, A4];

%% E and D
for i = 0:r
    E(i+1).E = computeE(s-i, i);
end

for i = 0:r
    D(i+1).D = computeD(s, i);
end

%% X CONTROL POINTS
EDM = [];
for i = 0:r-1
    EDM = [EDM; E(i+1).E*D(i+1).D];
end

%% Z CONTROL POINTS
At1 = kron(A1, eye(s+1));
At2 = kron(A2, eye(s+1));
Dz = computeD(s, 1);
Ez = computeE(s-1, 1);
Dt = [];
for i = 1:n-r
    Dt = blkdiag(Dt,Dz);
end

Et = [];
for i = 1:n-r
    Et = blkdiag(Et,Ez);
end

Maux = [eye((n-r)*s), -Dt; Et, -At1];
Maux1 = inv(Maux);
AtAux = Maux1(end-((n-r)*(s+1)-1):end, end-((n-r)*(s+1)-1):end)*At2*EDM;

for i = 1:n-r
    At3(i).A = AtAux( (i-1)*(s+1)+1:(i-1)*(s+1) + s+1, :);
end

%% U
Az = 0;
for i = 1:n-r
    Az = Az + a3(i)*At3(i).A;
end

Ax = 0;
for i = 0:r
    Ax = Ax + a4(i+1)*E(i+1).E * D(i+1).D;
end

U = (Ax - Az)*Y;

%% SIMULATION
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
z(:, 1) = z0;
for tt = 1:length(t)
    z(:, tt+1) = z(:, tt) + (A1*z(:, tt) + A2*x(:, tt))*dt;
    x(:, tt+1) = x(:, tt) + (A3*z(:, tt) + A4*x(:, tt) + B*u(tt))*dt;
end


plot(y)
hold on
plot(x(1, :))
hold off