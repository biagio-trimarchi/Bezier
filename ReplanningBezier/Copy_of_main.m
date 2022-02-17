% REPLANNING BEZIER
% Notes

%% Clear workspace
clear 
clc 
addpath('./Functions')

% generateBezierCompositionFun(3,3);


%% Parameters
n = 3;          % Spatial Curve Order
m = 3;          % Temporal Parametrization Order
d = 2;          % Space Order

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

Pd(:, 1) = [0.5; 0.0];
Pd(:, 2) = [0.2; 1.0];
Pd(:, 3) = [0.5; 2.4];
Pd(:, 4) = [1.2; 1.1];

% Temporal Control Points
Ud(1) = 0;
Ud(2) = 0.2;
Ud(3) = 0.7;
Ud(4) = 1;

tic
Cd = BezierComposition33(Pd(:, 1), Pd(:, 2), Pd(:, 3), Pd(:, 4), Ud(1), Ud(2), Ud(3), Ud(4))
toc

%% Bezier Curves
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

Cd = BezierComposition33(Pd(:, 1), Pd(:, 2), Pd(:, 3), Pd(:, 4), Ud(1), Ud(2), Ud(3), Ud(4));

%% Bezier Curves
bezierCurved2 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurved2(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Ud(k);
    end
    for k = 1:n+1
        bezierCurved2(:, tt) = bezierCurved2(:, tt) + bernsteinPol(n, k-1, u)*Pd(:, k);
    end
end
%% Generate Obstacle Path
% Spatial Control Points

% Trajectory 1
Po(:, 1) = [2.0; 0.0];
Po(:, 2) = [1.5; 1.0];
Po(:, 3) = [1.2; 2.4];
Po(:, 4) = [0.5; 1.1];

% Temporal Control Points
Uo(1) = 0;
Uo(2) = 0.1;
Uo(3) = 0.8;
Uo(4) = 1;

Co = BezierComposition33(Po(:, 1), Po(:, 2), Po(:, 3), Po(:, 4), Uo(1), Uo(2), Uo(3), Uo(4));

%% Bezier Curves
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

Co = BezierComposition33(Po(:, 1), Po(:, 2), Po(:, 3), Po(:, 4), Uo(1), Uo(2), Uo(3), Uo(4));

%% Bezier Curves
bezierCurveo2 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveo2(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Uo(k);
    end
    for k = 1:n+1
        bezierCurveo2(:, tt) = bezierCurveo2(:, tt) + bernsteinPol(n, k-1, u)*Po(:, k);
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
    Avu(k+1, k+1:k+d) = [-1, 1];
    Avd(k+1, k+1:k+d) = [1, -1];
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
bau = ones(n*m-1, 1)*aMax*(t2^2)/((n*m)*(n*m-1));
bad = ones(n*m-1, 1)*aMax*(t2^2)/((n*m)*(n*m-1));

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
x0 = [reshape(Pd, [], 1); Ud];
% x = [reshape(Po, [], 1); Uo];
% J = cost1(x, n, m, d, Pd, Ud)
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
opts = optimoptions('fmincon','Algorithm','sqp');

xNew = x0;
tic
[xNew, fval] = fmincon(J, x0, [], [], Aeq, beq, lb, ub, constraints, opts);
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
Cnew = BezierComposition33(PNew(:, 1), PNew(:, 2), PNew(:, 3), PNew(:, 4), UNew(1), UNew(2), UNew(3), UNew(4));

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
% Plot Difference and Norm Of Square Distance
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
% plot(PNew(1, :), PNew(2, :), 'o');
% hold on
% plot(bezierCurveNew(1, :), bezierCurveNew(2, :));
% hold off
% title('New Drone Trajectory')
% legend('Spatial Control Points', 'Trajectory')
% 
% Plot new norm of squared distance
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
% title('Norm of Difference Between Trajectorie% figure(9)
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
% legend('Trajectory')


%% Plot animated trajectories

% F(length(t)) = struct('cdata',[],'colormap',[]);
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
% 
% for tt =1:length(t)
%     figure(12)
%     plot(bezierCurveo(1, :), bezierCurveo(2, :));
%     hold on
%     plot(bezierCurved(1, :), bezierCurved(2, :));
%     plot(bezierCurveNew(1, :), bezierCurveNew(2, :));
%     
%      
%     %plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
%     %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
%     %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
%     plotCircle(bezierCurveNew(1, tt), bezierCurveNew(2, tt), safeDist)
%     plotCircle(bezierCurveNew(1, tt), bezierCurveNew(2, tt), rDetect)
%     plot(bezierCurveNew(1, tt), bezierCurveNew(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
%     plot(bezierCurveo2(1, tt), bezierCurveo2(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
%     hold off 
%     
%     xlim([-1, 4])
%     ylim([-1, 4])
%     title("Simulation New Trajectory")
%     legend("Obstacle Trajectory", "Old Trjajectory", "New Trajectory", "Safe Circle", "Detect Distance", 'Location', 'southeast')
%     drawnow 
%     F(tt) = getframe;
% end
