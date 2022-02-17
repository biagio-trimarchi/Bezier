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

rDetect = 1;        % Detection Radius
safeDist = 0.25;    % Safe Distance 
vMax = 10;          % Maximum Velocity
aMax = 5;           % Maximum Acceleration
tF = 50;            % Total Time
t1 = 17;            % Detection Time
t2 = tF-t1;         % Remaining Time
t = 0:0.01:1;       % Time Vector

%% Generate Drone Path
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

%% Decompose Path
[Ud1, Ud2] = bezierDivision(Ud(1, :), t1/tF, m);
Ud2 = (Ud2-Ud2(1))/(1-Ud2(1));

uDiv = 0;
for k = 1:m+1
    uDiv = uDiv + bernsteinPol(m, k-1, t1/tF)*Ud(k);
end
uDiv1 = uDiv;

[Px1, Px2] = bezierDivision(Pd(1, :), uDiv, n);
[Py1, Py2] = bezierDivision(Pd(2, :), uDiv, n);

Ud = Ud2;
Pd = [Px2'; Py2'];

Cd = BezierComposition(Pd(:, 1), Pd(:, 2), Pd(:, 3), Pd(:, 4), Pd(:, 5), Pd(:, 6), Ud(1), Ud(2), Ud(3), Ud(4));
Cd = reshape(Cd, 2, []);

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

%% Decompose Path
[Uo1, Uo2] = bezierDivision(Uo(1, :), t1/tF, m);
Uo2 = (Uo2-Uo2(1))/(1-Uo2(1));

uDiv = 0;
for k = 1:m+1
    uDiv = uDiv + bernsteinPol(m, k-1, t1/tF)*Uo(k);
end
[Px1, Px2] = bezierDivision(Po(1, :), uDiv, n);
[Py1, Py2] = bezierDivision(Po(2, :), uDiv, n);

Uo = Uo2;
Po = [Px2'; Py2'];

Co = BezierComposition(Po(:, 1), Po(:, 2), Po(:, 3), Po(:, 4), Po(:, 5), Po(:, 6), ...
                       Uo(1), Uo(2), Uo(3), Uo(4));
Co = reshape(Co, 2, []);

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
% Set Initial Condition
x0 = reshape(Cd, 1, []);
%J = cost2(x0, n*m, d, Cd);
%[c, ceq] = nonLinearConstraints2(x0, n*m, d, Co, Q, safeDist);

% Inequality Constraints
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

% Equality Constraints
Aeq(1:d, :) = [eye(d), zeros(d, (n*m)*d)];      % Initial Position
Aeq(d+1:d+d, :) = [zeros(d, (n*m)*d), eye(d)];  % Final Position
Aeq(d+d+1:d+d+d, :) = [-eye(d), eye(d), zeros(d, (n*m-1)*d)];   % Initial Velocity
Aeq(d+d+1:d+d+d, :) = [zeros(d, (n*m-1)*d), -eye(d), eye(d)];   % Final Velocity

beq(1:d, 1) = Cd(:, 1);
beq(d+1:d+d, 1) = Cd(:, end);
beq(d+d+1:d+d+d, 1) = Cd(:, 2) - Cd(:, 1);
beq(d+d+1:d+d+d, 1) = Cd(:, end) - Cd(:, end-1);

% Function Handle for cost and non linear constraints
aux = reshape(Cd, 1, []);
J = @(x) cost2(x, n*m, d, aux);
constraints = @(x) nonLinearConstraints2Old(x, n*m, d, Co, Q, safeDist);
options = optimoptions('fmincon','SpecifyObjectiveGradient',true, 'Algorithm','sqp');

tic
[xNew, fval] = fmincon(J, x0, A, b, Aeq, beq, [], [], constraints);
toc
%% New Drone Trajectory
CNew = reshape(xNew, d, []);

% Bezier Curve
bezierCurveNew = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveNew(:, tt) = zeros(d, 1);
    for k = 1:(n*m)+1
        bezierCurveNew(:, tt) = bezierCurveNew(:, tt) + bernsteinPol(n*m, k-1, t(tt))*CNew(:, k);
    end
end

CdiffNew = CNew - Co;

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
% % Plot Drone Trajectory
% figure(1)
% plot(Pd(1, :), Pd(2, :), 'o');
% hold on
% plot(bezierCurved(1, :), bezierCurved(2, :));
% hold off
% title('Drone Initial Trajectory')
% legend('Spatial Control Points', 'Trajectory')
% 
% figure(2)
% plot(Cd(1, :), Cd(2, :), 'x');
% hold on
% plot(bezierCurve2d(1, :), bezierCurve2d(2, :));
% hold off
% title('Drone Initial Trajectory')
% legend('Expanded Control Points', 'Trajectory')
% 
% % Plot Obstacle Trajectory
% figure(3)
% plot(Po(1, :), Po(2, :), 'o');
% hold on
% plot(bezierCurveo(1, :), bezierCurveo(2, :));
% hold off
% title('Obstacle Trajectory')
% legend('Spatial Control Points', 'Trajectory')
% 
% figure(4)
% plot(Co(1, :), Co(2, :), 'x');
% hold on
% plot(bezierCurve2o(1, :), bezierCurve2o(2, :));
% hold off
% title('Obstacle Trajectory')
% legend('Expanded Control Points', 'Trajectory')
% 
% % Plot Difference and Norm Of Square Distance
% figure(5)
% plot(Cdiff(1, :), Cdiff(2, :), 'x');
% hold on
% plot(bezierCurveDiff(1, :), bezierCurveDiff(2, :));
% hold off
% title('Difference Between Trajectories')
% legend('Control Points', 'Trajectory')
% 
% figure(6)
% plot(CdiffNorm, 'x');
% hold on 
% plot(safeDist^2*ones(length(CdiffNorm), 1))
% hold off
% title('Norm of Difference Between Trajectories')
% legend('Control Points')
% 
% figure(7)
% plot(t, bezierCurveDiffNorm);
% hold on 
% plot(t, safeDist^2*ones(length(t), 1))
% hold off
% title('Norm of Difference Between Trajectories')
% legend('Trajectory')
% 
% % Plot New Drone Trajectory
% figure(8)
% plot(CNew(1, :), CNew(2, :), 'o');
% hold on
% plot(bezierCurveNew(1, :), bezierCurveNew(2, :));
% hold off
% title('New Drone Trajectory')
% legend('Spatial Control Points', 'Trajectory')
% 
% % Plot new norm of squared distance
% figure(9)
% plot(CdiffNormNew, 'x');
% hold on 
% plot(safeDist^2*ones(length(CdiffNormNew), 1))
% hold off
% title('Norm of Difference Between Trajectories')
% legend('Control Points')
% 
% figure(10)
% plot(t, bezierCurveDiffNormNew);
% hold on 
% plot(t, safeDist^2*ones(length(t), 1))
% hold off
% title('Norm of Difference Between Trajectories')
% legend('Trajectory')
% 
% 
% %% Plot animated trajectories
% 
% for tt =1:length(t)
%     figure(11)
%     plot(bezierCurveo(1, :), bezierCurveo(2, :));
%     hold on
%     plot(bezierCurved(1, :), bezierCurved(2, :));
%     plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
%     plot(bezierCurveo(1, tt), bezierCurveo(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
%     plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
%     hold off 
%     
%     xlim([-1, 4])
%     ylim([-1, 4])
%     title("Simulation Original Trajectory")
%     drawnow 
% end
% 
% for tt =1:length(t)
%     figure(12)
%     plot(bezierCurveo(1, :), bezierCurveo(2, :));
%     hold on
%     plot(bezierCurveNew(1, :), bezierCurveNew(2, :));
%     plot(bezierCurved(1, :), bezierCurved(2, :));
%     plot(bezierCurveNew(1, tt), bezierCurveNew(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
%     plot(bezierCurveo(1, tt), bezierCurveo(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
%     plotCircle(bezierCurveNew(1, tt), bezierCurveNew(2, tt), safeDist)
%     hold off 
%     
%     xlim([-1, 4])
%     ylim([-1, 4])
%     title("Simulation New Trajectory")
%     drawnow 
% end
