% REPLANNING BEZIER
% Notes

%% Clear workspace
clear 
clc 
addpath('./Functions')

%% PARAMETERS
n = 5;
m = 3;
d = 2;
dT = 2;
N = 2;
Tf = N*dT;
t = 0:0.01:1;   % Time Vector
dSafe = 1;

%% COMPUTE Q
Q = computeQ(n*m);

%% DRONES CONTROL POINTS

P1 = [0,0,0,1,2,3; 
     0,0,0,3,1,4;
     3,4,5,  6,7,7;
     4,7,15, 7,4,4];
 
U1 = [0, 0.4, 0.7, 1;
      0, 0.3, 0.5, 1];

P2 = [3, 3, 3, 2, 1, 3.5;  
      0, 0, 0, 1, 2, 4.5;
      3.5, 5.5, 11, 3, 2, 4;
      4.5, 6.5, 7, 2, 1, 0];
  
U2 = [0, 0.4, 0.7, 1;
      0, 0.3, 0.5, 1];
  
%% OPTIMIZATION PROBLEM

% Equality
AP0 = [1, zeros(1, n)];
APF = [zeros(1, n), 1];
AU0 = [1, zeros(1, m)];
AUF = [zeros(1, m), 1];

Aeq = blkdiag(AP0, APF, AP0, APF);
Aeq = [Aeq; blkdiag([APF, -AP0], [APF, -AP0])];
Aeq = blkdiag(Aeq, [AU0; AUF], [AU0; AUF]);

beq = [P1(1, 1);
       P1(3, end);
       P1(2, 1);
       P1(4, end);
       0;
       0;
       0;
       1;
       0;
       1;];

%Dynamic

% Bounds

lb = [ones(2*d*(n+1), 1)*(-inf); zeros(2*(m+1), 1)];
ub = [ones(2*d*(n+1), 1)*(inf); ones(2*(m+1), 1)];
% Parameters
Co = BezierComposition(P2(1:2, 1),P2(1:2, 2),P2(1:2, 3),P2(1:2, 4),P2(1:2, 5),P2(1:2 ,6),U2(:, 1),U2(:, 2),U2(:, 3),U2(: ,4));
Co = reshape(Co, d, []);

%Functions
J = @(x) costPU(x, n ,d, P1(1:2, n+1));
constr = @(x) nlCnstPU(x, Co, n, m, d, dSafe, Q, P1, U1);


%Optimization
x0 = [P1(1, :), P1(3, :), P1(2, :), P1(4, :), U1(1, :), U2(2, :)];
nlCnstPU(x0, Co, n, m, d, dSafe, Q, P1, U1);

opts = optimoptions('fmincon','Algorithm','sqp');
tic
[x, fVal] = fmincon(J, x0, [], [], Aeq, beq, lb, ub, constr, opts);
toc

P3 = [x(1: n+1); x(2*(n+1)+1:3*(n+1)); x(n+2: 2*(n+1)); x(3*(n+1)+1: 4*(n+1))];
U3 = [x(4*(n+1)+1:4*(n+1)+m+1);
      x(4*(n+1)+m+2:end)];

%% CURVES

bezierCurve11 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve11(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U1(1, k);
    end
    for k = 1:n+1
        bezierCurve11(:, tt) = bezierCurve11(:, tt) + bernsteinPol(n, k-1, u)*P1(1:2, k);
    end
end

bezierCurve12 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve12(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U1(2, k);
    end
    for k = 1:n+1
        bezierCurve12(:, tt) = bezierCurve12(:, tt) + bernsteinPol(n, k-1, u)*P1(3:4, k);
    end
end

bezierCurve21 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve21(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U2(1, k);
    end
    for k = 1:n+1
        bezierCurve21(:, tt) = bezierCurve21(:, tt) + bernsteinPol(n, k-1, u)*P2(1:2, k);
    end
end

bezierCurve22 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve22(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U2(2, k);
    end
    for k = 1:n+1
        bezierCurve22(:, tt) = bezierCurve22(:, tt) + bernsteinPol(n, k-1, u)*P2(3:4, k);
    end
end

bezierCurve31 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve31(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U3(1, k);
    end
    for k = 1:n+1
        bezierCurve31(:, tt) = bezierCurve31(:, tt) + bernsteinPol(n, k-1, u)*P3(1:2, k);
    end
end

bezierCurve32 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve32(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U2(2, k);
    end
    for k = 1:n+1
        bezierCurve32(:, tt) = bezierCurve32(:, tt) + bernsteinPol(n, k-1, u)*P3(3:4, k);
    end
end

%% PLOTS
figure(1)
plot(bezierCurve11(1, :), bezierCurve11(2, :));
hold on
plot(bezierCurve12(1, :), bezierCurve12(2, :));
plot(bezierCurve21(1, :), bezierCurve21(2, :));
plot(bezierCurve22(1, :), bezierCurve22(2, :));
plot(P2(1, :), P2(2, :), 'o');
hold off

figure(2)
plot(bezierCurve31(1, :), bezierCurve31(2, :));
hold on
plot(bezierCurve32(1, :), bezierCurve32(2, :));
plot(bezierCurve21(1, :), bezierCurve21(2, :));
plot(bezierCurve22(1, :), bezierCurve22(2, :));
plot(P2(1, :), P2(2, :), 'o');
hold off

for tt =1:length(t)
    figure(12)
    plot(bezierCurve11(1, :), bezierCurve11(2, :));
    hold on
    plot(bezierCurve21(1, :), bezierCurve21(2, :));
    
     
    %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
    plotCircle(bezierCurve11(1, tt), bezierCurve11(2, tt), dSafe)
    plot(bezierCurve11(1, tt), bezierCurve11(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    plot(bezierCurve21(1, tt), bezierCurve21(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
    hold off 
    
    xlim([-1, 11])
    ylim([-1, 11])
    title("Simulation New Trajectory")
    drawnow 
end

for tt =1:length(t)
    figure(12)
    plot(bezierCurve31(1, :), bezierCurve31(2, :));
    hold on
    plot(bezierCurve21(1, :), bezierCurve21(2, :));
    
     
    %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
    plotCircle(bezierCurve31(1, tt), bezierCurve31(2, tt), dSafe)
    plot(bezierCurve31(1, tt), bezierCurve31(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    plot(bezierCurve21(1, tt), bezierCurve21(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
    hold off 
    
    xlim([-1, 11])
    ylim([-1, 11])
    title("Simulation New Trajectory")
    drawnow 
end