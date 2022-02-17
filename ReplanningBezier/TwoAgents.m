% REPLANNING BEZIER
% Notes

%% Clear workspace
clear 
clc 
addpath('./Functions')

%generateBezierCompositionFun(5,3);


%% Parameters
n = 5;          % Spatial Curve Order
m = 3;          % Temporal Parametrization Order
d = 2;          % Space Order
t = 0:0.01:1;   % Time Vector

safeDist = 0.25;    % Safe Distance 
vMax = 10;          % Maximum Velocity
aMax = 5;          % Maximum Acceleration
tF = 10;            % Total Time

%% Generate Drone1 Path
% Spatial Control Points

Pd(:, 1) = [0.0; 0.0];
Pd(:, 2) = [0.2; 1.0];
Pd(:, 3) = [0.5; 2.4];
Pd(:, 4) = [1.2; 1.1];
Pd(:, 5) = [1.5; 2.0];
Pd(:, 6) = [2.0; 3.0];

% Temporal Control Points
Ud(1) = 0;
Ud(2) = 0.2;
Ud(3) = 0.7;
Ud(4) = 1;

% %% Decompose Path
% [Ud1, Ud2] = bezierDivision(Ud(1, :), t1/tF, m);
% Ud2 = (Ud2-Ud2(1))/(1-Ud2(1));
% 
% uDiv = 0;
% for k = 1:m+1
%     uDiv = uDiv + bernsteinPol(m, k-1, t1/tF)*Ud(k);
% end
% uDiv1 = uDiv;
% [Px1, Px2] = bezierDivision(Pd(1, :), uDiv, n);
% [Py1, Py2] = bezierDivision(Pd(2, :), uDiv, n);
% 
% Ud = Ud2;
% Pd = [Px2'; Py2'];

%% Generate Drone2 Path
% Spatial Control Points

% Trajectory 1
Po(:, 1) = [2.0; 0.0];
Po(:, 2) = [1.5; 1.0];
Po(:, 3) = [1.2; 2.4];
Po(:, 4) = [0.5; 1.1];
Po(:, 5) = [0.2; 2.0];
Po(:, 6) = [0.0; 3.0];

% Trajectory 2
% Po(:, 1) = [2.0; 0.0];
% Po(:, 2) = [1.5; 0.5];
% Po(:, 3) = [1.2; 0.2];
% Po(:, 4) = [0.5; 0.4];
% Po(:, 5) = [0.2; 0.4];
% Po(:, 6) = [0.0; 0.0];

% Temporal Control Points
Uo(1) = 0;
Uo(2) = 0.1;
Uo(3) = 0.8;
Uo(4) = 1;

%% Compute Matrix for Square Norm
Q = computeQ(n*m);

%% Solve optimization Problem
% Roba
% Velocità
Avu = zeros(n*m, n*m+1);
Avd = zeros(n*m, n*m+1);
bvu = ones(n*m, 1)*vMax*tF/(n*m);
bvd = ones(n*m, 1)*vMax*tF/(n*m);

for k = 0:n*m-1
    Avu(k+1, k+1:k+2) = [-1, 1];
    Avd(k+1, k+1:k+2) = [1, -1];
end
Avu = kron(Avu, eye(d));
bvu = kron(bvu, ones(d, 1));
Avd = kron(Avd, eye(d));
bvd = kron(bvd, ones(d, 1));

Av = [Avu; Avd];
bv = [bvu; bvd];

% Accelerazione
Aau = zeros(n*m-1, n*m+1);
Aad = zeros(n*m-1, n*m+1);
bau = ones(n*m-1, 1)*aMax*(tF^2)/((n*m)*(n*m-1));
bad = ones(n*m-1, 1)*aMax*(tF^2)/((n*m)*(n*m-1));

for k = 0:n*m-2
    Aau(k+1, k+1:k+3) = [1, -2, 1];
    Aad(k+1, k+1:k+3) = [-1, 2, -1];
end
Aau = kron(Aau, eye(d));
bau = kron(bau, ones(d, 1));
Aad = kron(Aad, eye(d));
bad = kron(bad, ones(d, 1));

Aa = [Aau; Aad];
ba = [bau; bad];

A = [Av; Aa];
b = [bv; ba];

% Set Initial Condition
x0 = [reshape(Pd, [], 1); Ud'; reshape(Po, [], 1); Uo'];
J = cost2Agents(x0, n, m, d, Pd, Ud, Po, Uo);
[c, ceq] = nonLinearConstraints2Agents(x0, n, m, d, Q, Pd, Ud, Po, Uo, safeDist, A, b);

% Linear Constraints
Aeq(1:d, :) = [eye(d), zeros(d, n*d+m+1)];
Aeq(d+1:d+d, :) = [zeros(d, n*d), eye(d), zeros(d, m+1)];
Aeq(d+d+1, :) = [zeros(1, (n+1)*d), 1, zeros(1, m)];
Aeq(d+d+2, :) = [zeros(1, (n+1)*d+m), 1];

Aeq = blkdiag(Aeq, Aeq);

beq(1:d, 1) = Pd(:, 1);
beq(d+1:d+d, 1) = Pd(:, end);
beq(d+d+1, 1) = 0;
beq(d+d+2, 1) = 1;
beq(d+d+3:d+d+2+d, 1) = Po(:, 1);
beq(d+d+2+d+1:d+d+2+d+d, 1) = Po(:, end);
beq(d+d+2+d+d+1, 1) = 0;
beq(d+d+2+d+d+2, 1) = 1;

% Bounds
lb = [-inf*ones((n+1)*d, 1); zeros(m+1, 1); -inf*ones((n+1)*d, 1); zeros(m+1, 1)];
ub = [inf*ones((n+1)*d, 1); ones(m+1, 1); inf*ones((n+1)*d, 1); ones(m+1, 1)];

% Function Handle for cost and non linear constraints
J = @(x) cost2Agents(x, n, m, d, Pd, Ud, Po, Uo);
constraints = @(x) nonLinearConstraints2Agents(x, n, m, d, Q, Pd, Ud, Po, Uo, safeDist, A, b);

opts = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',10000);

xNew = x0;
tic
[xNew, fval] = fmincon(J, x0, [], [], Aeq, beq, lb, ub, constraints, opts);
toc

P1 = xNew(1:(n+1)*d);
P1 = reshape(P1, d, []);
U1 = xNew((n+1)*d+1:(n+1)*d+m+1);

P2 = xNew((n+1)*d+m+2:(n+1)*d+m +(n+1)*d+1);
P2 = reshape(P2, d, []);
U2 = xNew((n+1)*d+m +(n+1)*d+2:end);

bezierCurve1 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve1(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U1(k);
    end
    for k = 1:n+1
        bezierCurve1(:, tt) = bezierCurve1(:, tt) + bernsteinPol(n, k-1, u)*P1(:, k);
    end
end

bezierCurve2 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve2(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*U2(k);
    end
    for k = 1:n+1
        bezierCurve2(:, tt) = bezierCurve2(:, tt) + bernsteinPol(n, k-1, u)*P2(:, k);
    end
end

bezierCurve3 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve3(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Ud(k);
    end
    for k = 1:n+1
        bezierCurve3(:, tt) = bezierCurve3(:, tt) + bernsteinPol(n, k-1, u)*Pd(:, k);
    end
end

bezierCurve4 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve4(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Uo(k);
    end
    for k = 1:n+1
        bezierCurve4(:, tt) = bezierCurve4(:, tt) + bernsteinPol(n, k-1, u)*Po(:, k);
    end
end

%% Plot animated trajectories
% 
% for tt =1:length(t)
%     figure(11)
%     plot(bezierCurveo(1, :), bezierCurveo(2, :));
%     hold on
%     plot(bezierCurved(1, :), bezierCurved(2, :));
%     plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
%     plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
%     plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
%     plot(bezierCurveo(1, tt), bezierCurveo(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
%     hold off 
%     
%     xlim([-1, 4])
%     ylim([-1, 4])
%     title("Simulation Original Trajectory")
%     %text(0.5,0.5,sprintf('%f',(tt-1)*tF/length(t)))
%     legend("Obstacle Trajectory", "Drone Trajectory", "Safe Circle", "Detect Distance", 'Location', 'southeast')
%     drawnow 
%     
% end

for tt =1:length(t)
    figure(11)
    plot(bezierCurve3(1, :), bezierCurve3(2, :));
    hold on
    plot(bezierCurve4(1, :), bezierCurve4(2, :));
    
     
    %plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
    %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
    plotCircle(bezierCurve3(1, tt), bezierCurve3(2, tt), safeDist)
    plot(bezierCurve3(1, tt), bezierCurve3(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    plot(bezierCurve4(1, tt), bezierCurve4(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
    hold off 
    
    xlim([-1, 4])
    ylim([-1, 4])
    title("Simulation New Trajectory")
    drawnow 
end

for tt =1:length(t)
    figure(12)
    plot(bezierCurve1(1, :), bezierCurve1(2, :));
    hold on
    plot(bezierCurve2(1, :), bezierCurve2(2, :));
    
     
    %plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
    %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
    plotCircle(bezierCurve1(1, tt), bezierCurve1(2, tt), safeDist)
    plot(bezierCurve1(1, tt), bezierCurve1(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    plot(bezierCurve2(1, tt), bezierCurve2(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
    hold off 
    
    xlim([-1, 4])
    ylim([-1, 4])
    title("Simulation New Trajectory")
    drawnow 
end