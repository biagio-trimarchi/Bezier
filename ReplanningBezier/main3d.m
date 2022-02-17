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
d = 3;          % Space Order
t = 0:0.01:1;   % Time Vector

safeDist = 0.25;        % Safe Distance 
vMax = 20;              % Maximum Velocity
aMax = 10;               % Maximum Acceleration
tF = 20;                % Total Time

%% Generate Drone Path
% Spatial Control Points

Pd(:, 1) = [0.0; 0.0; 0.0];
Pd(:, 2) = [0.2; 1.0; 0.0];
Pd(:, 3) = [0.5; 2.4; 0.0];
Pd(:, 4) = [1.2; 1.1; 0.0];
Pd(:, 5) = [1.5; 2.0; 0.0];
Pd(:, 6) = [2.0; 3.0; 1.0];

% Temporal Control Points
Ud(1) = 0;
Ud(2) = 0.2;
Ud(3) = 0.7;
Ud(4) = 1;

Cd = BezierComposition(Pd(:, 1), Pd(:, 2), Pd(:, 3), Pd(:, 4), Pd(:, 5), ...
                      Pd(:, 6), Ud(1), Ud(2), Ud(3), Ud(4));
Cd = reshape(Cd, d, []);

% Bezier Curves
bezierCurved = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurved(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Ud(k);
    end
    for k = 1:n+1
        bezierCurved(:, tt) = bezierCurved(:, tt) + bernsteinPol(n, k-1, u)*Pd(:, k);
    end
end


bezierCurve2d = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve2d(:, tt) = zeros(d, 1);
    for k = 1:(m*n)+1
        bezierCurve2d(:, tt) = bezierCurve2d(:, tt) + bernsteinPol(m*n, k-1, t(tt))*Cd(:, k);
    end
end

%% Generate Obstacle Path
% Spatial Control Points

% Trajectory 1
Po(:, 1) = [2.0; 0.0; 0.0];
Po(:, 2) = [1.5; 1.0; 0.0];
Po(:, 3) = [1.2; 2.4; 0.0];
Po(:, 4) = [0.5; 1.1; 0.0];
Po(:, 5) = [0.2; 2.0; 0.0];
Po(:, 6) = [0.0; 3.0; 1.0];

% Trajectory 2
% Po(:, 1) = [2.0; 0.0; 0.0];
% Po(:, 2) = [1.5; 0.5; 0.0];
% Po(:, 3) = [1.2; 0.2; 0.0];
% Po(:, 4) = [0.5; 0.4; 0.0];
% Po(:, 5) = [0.2; 0.4; 0.0];
% Po(:, 6) = [0.0; 0.0; 0.0];

% Temporal Control Points
Uo(1) = 0;
Uo(2) = 0.1;
Uo(3) = 0.8;
Uo(4) = 1;

Co = BezierComposition(Po(:, 1), Po(:, 2), Po(:, 3), Po(:, 4), Po(:, 5), ...
                      Po(:, 6), Uo(1), Uo(2), Uo(3), Uo(4));
Co = reshape(Co, d, []);

% Bezier Curves
bezierCurveo = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveo(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Uo(k);
    end
    for k = 1:n+1
        bezierCurveo(:, tt) = bezierCurveo(:, tt) + bernsteinPol(n, k-1, u)*Po(:, k);
    end
end


bezierCurve2o = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve2o(:, tt) = zeros(d, 1);
    for k = 1:(m*n)+1
        bezierCurve2o(:, tt) = bezierCurve2o(:, tt) + bernsteinPol(m*n, k-1, t(tt))*Co(:, k);
    end
end

%% Compute Matrix for Square Norm
Q = computeQ(n*m);

%% Compute Square Norm
Cdiff = Cd - Co;

CdiffNorm = zeros(1,2*(n*m)+1);
for k = 0:2*n*m
    CdiffNorm(k+1) = 0;
    for j = 1:d
    CdiffNorm(k+1) = CdiffNorm(k+1) + ...
                     Cdiff(j, :)*Q(:, :, k+1)*Cdiff(j, :)';

    end
end

bezierCurveDiff = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveDiff(:, tt) = zeros(d, 1);
    for k = 1:(m*n)+1
        bezierCurveDiff(:, tt) = bezierCurveDiff(:, tt) + bernsteinPol(m*n, k-1, t(tt))*Cdiff(:, k);
    end
end

bezierCurveDiffNorm = zeros(1, length(t));
for tt = 1:length(t)
    bezierCurveDiffNorm(tt) = 0;
    for k = 1:(2*m*n)+1
        bezierCurveDiffNorm(tt) = bezierCurveDiffNorm(tt) + bernsteinPol(2*m*n, k-1, t(tt))*CdiffNorm(:, k);
    end
end

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
x0 = [reshape(Pd, [], 1); Ud'];
%J = cost1(x, n, m, d, Pd, Ud);
%[c, ceq] = nonLinearConstraints(x, n, m, d, Q, Co);

% Linear Constraints
Aeq(1:d, :) = [eye(d), zeros(d, n*d+m+1)];
Aeq(d+1:d+d, :) = [zeros(d, n*d), eye(d), zeros(d, m+1)];
Aeq(d+d+1, :) = [zeros(1, (n+1)*d), 1, zeros(1, m)];
Aeq(d+d+2, :) = [zeros(1, (n+1)*d+m), 1];

beq(1:d, 1) = Pd(:, 1);
beq(d+1:d+d, 1) = Pd(:, end);
beq(d+d+1, 1) = 0;
beq(d+d+2, 1) = 1;

% Bounds
lb = [-inf*ones((n+1)*d, 1); zeros(m+1, 1)];
ub = [inf*ones((n+1)*d, 1); ones(m+1, 1)];

% Function Handle for cost and non linear constraints
J = @(x) cost1(x, n, m, d, Pd, Ud);
constraints = @(x) nonLinearConstraints(x, n, m, d, Q, Co, Pd, Ud, safeDist, A, b);

xNew = x0;
tic
[xNew, fval] = fmincon(J, x0, [], [], Aeq, beq, lb, ub, constraints);
toc
%% New Drone Trajectory
PNew = xNew(1:(n+1)*d);
PNew = reshape(PNew, d, []);
UNew = xNew((n+1)*d+1:end);

% Bezier Curve
bezierCurveNew = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveNew(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*UNew(k);
    end
    for k = 1:n+1
        bezierCurveNew(:, tt) = bezierCurveNew(:, tt) + bernsteinPol(n, k-1, u)*PNew(:, k);
    end
end

%% New Norm of Squared Distance 
Cnew = BezierComposition(PNew(:, 1), PNew(:, 2), PNew(:, 3), PNew(:, 4), PNew(:, 5), ...
                      PNew(:, 6), UNew(1), UNew(2), UNew(3), UNew(4));
Cnew = reshape(Cnew, d, []);

CdiffNew = Cnew - Co;

CdiffNormNew = zeros(1,2*(n*m)+1);
for k = 0:2*n*m
    CdiffNormNew(k+1) = 0;
    for j = 1:d
    CdiffNormNew(k+1) = CdiffNormNew(k+1) + ...
                     CdiffNew(j, :)*Q(:, :, k+1)*CdiffNew(j, :)';
    end
end

bezierCurveDiffNormNew = zeros(1, length(t));
for tt = 1:length(t)
    bezierCurveDiffNormNew(tt) = 0;
    for k = 1:(2*m*n)+1
        bezierCurveDiffNormNew(tt) = bezierCurveDiffNormNew(tt) + bernsteinPol(2*m*n, k-1, t(tt))*CdiffNormNew(:, k);
    end
end

%% Plots
% Plot Drone Trajectory
figure(1)
plot3(Pd(1, :), Pd(2, :), Pd(3, :), 'o');
hold on
plot3(bezierCurved(1, :), bezierCurved(2, :), bezierCurved(2, :));
hold off
title('Drone Initial Trajectory')
legend('Spatial Control Points', 'Trajectory')

figure(2)
plot3(Cd(1, :), Cd(2, :), Cd(3, :), 'x');
hold on
plot3(bezierCurve2d(1, :), bezierCurve2d(2, :), bezierCurve2d(3, :));
hold off
title('Drone Initial Trajectory')
legend('Expanded Control Points', 'Trajectory')

% Plot Obstacle Trajectory
figure(3)
plot3(Po(1, :), Po(2, :), Po(3, :), 'o');
hold on
plot3(bezierCurveo(1, :), bezierCurveo(2, :), bezierCurveo(3, :));
hold off
title('Obstacle Trajectory')
legend('Spatial Control Points', 'Trajectory')

figure(4)
plot3(Co(1, :), Co(2, :), Co(3, :), 'x');
hold on
plot3(bezierCurve2o(1, :), bezierCurve2o(2, :), bezierCurve2o(3, :));
hold off
title('Obstacle Trajectory')
legend('Expanded Control Points', 'Trajectory')

% Plot Difference and Norm Of Square Distance
figure(5)
plot3(Cdiff(1, :), Cdiff(2, :), Cdiff(3, :), 'x');
hold on
plot3(bezierCurveDiff(1, :), bezierCurveDiff(2, :), bezierCurveDiff(3, :));
hold off
title('Difference Between Trajectories')
legend('Control Points', 'Trajectory')

figure(6)
plot(CdiffNorm, 'x');
hold on 
plot(safeDist^2*ones(length(CdiffNorm), 1))
hold off
title('Norm of Difference Between Trajectories')
legend('Control Points')

figure(7)
plot(t, bezierCurveDiffNorm);
hold on 
plot(t, safeDist^2*ones(length(t), 1))
hold off
title('Norm of Difference Between Trajectories')
legend('Trajectory')

% Plot New Drone Trajectory
figure(8)
plot3(PNew(1, :), PNew(2, :), PNew(3, :), 'o');
hold on
plot3(bezierCurveNew(1, :), bezierCurveNew(2, :), bezierCurveNew(3, :));
hold off
title('New Drone Trajectory')
legend('Spatial Control Points', 'Trajectory')

% Plot new norm of squared distance
figure(9)
plot(CdiffNormNew, 'x');
hold on 
plot(safeDist^2*ones(length(CdiffNormNew), 1))
hold off
title('Norm of Difference Between Trajectories')
legend('Control Points')

figure(10)
plot(t, bezierCurveDiffNormNew);
hold on 
plot(t, safeDist^2*ones(length(t), 1))
hold off
title('Norm of Difference Between Trajectories')
legend('Trajectory')


%% Plot animated trajectories

for tt =1:length(t)
    figure(11)
    plot3(bezierCurveo(1, :), bezierCurveo(2, :), bezierCurveo(3, :));
    hold on
    plot3(bezierCurved(1, :), bezierCurved(2, :), bezierCurved(3, :));
    plot3(bezierCurved(1, tt), bezierCurved(2, tt), bezierCurved(3, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    plot3(bezierCurveo(1, tt), bezierCurveo(2, tt), bezierCurveo(3, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
    plotSphere(bezierCurved(1, tt), bezierCurved(2, tt), bezierCurved(3, tt), safeDist)
    hold off 
    
    xlim([-1, 4])
    ylim([-1, 4])
    zlim([-2, 2])
    title("Simulation Original Trajectory")
    drawnow limitrate
end

for tt =1:length(t)
    figure(12)
    plot3(bezierCurveo(1, :), bezierCurveo(2, :), bezierCurveo(3, :));
    hold on
    plot3(bezierCurveNew(1, :), bezierCurveNew(2, :), bezierCurveNew(3, :));
    plot3(bezierCurveNew(1, tt), bezierCurveNew(2, tt), bezierCurveNew(3, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
    plot3(bezierCurveo(1, tt), bezierCurveo(2, tt), bezierCurveo(3, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
    plotSphere(bezierCurveNew(1, tt), bezierCurveNew(2, tt), bezierCurveNew(3, tt), safeDist)
    hold off 
    
    xlim([-1, 4])
    ylim([-1, 4])
    zlim([-1, 4])
    title("Simulation New Trajectory")
    drawnow limitrate
end
