%% %% %% AVOIDANCE CLUTTERED %% %% %%
% Notes

%% SETUP
clear 
clc 
addpath('./Functions')

%% PARAMETERS
n = 5;              % Spatial Beziér Curve Order
m = 3;              % Temporal Beziér Curve Order
d = 3;              % Space Dimension
r = 2;              % Shape Parameter
segments_num = 3;   % Number of Segments
t = 0:0.01:1;       % Time Vector (for visualization)

%% %% INITIAL MISSION PLANNING
%% MISSION PARAMETERS
start_position1 = [-5; 7; 7];                 % Start Position
target_position1 = [15; 6; 2];               % Target Position
start_position2 = [-5; 2; 7];                 % Start Position
target_position2 = [15; 2; 3];               % Target Position
segment_time = 5;                           % Time Allocated to each Segment
total_time = segments_num*segment_time;     % Total Mission Time

%% VEHICLE PARAMETERS
velocity_max = 100;         % Maximum Velocity
acceleration_max = 20;     % Maximum Acceleration
safeDist = 1.2;              % SafeDistance

%% SAFE REGIONS
% Safe Region 1
AE1x = eye(n+1);
AE1y = [-eye(n+1); eye(n+1)];
AE1z = [-eye(n+1); eye(n+1)];
bE1x = 5*ones(n+1, 1) - safeDist;
bE1y = [zeros(n+1, 1) - safeDist; 9*ones(n+1, 1) - safeDist];
bE1z = [zeros(n+1, 1) - safeDist; 9*ones(n+1, 1) - safeDist];

% Safe Region 2
AE2x = [-eye(n+1); eye(n+1)];
AE2y = [-eye(n+1); eye(n+1)];
AE2z = [-eye(n+1); eye(n+1)];
bE2x = [10*ones(n+1, 1) - safeDist; 20*ones(n+1, 1) - safeDist];
bE2y = [-3*ones(n+1, 1) - safeDist; 6*ones(n+1, 1) - safeDist];
bE2z = [-3*ones(n+1, 1) - safeDist; 6*ones(n+1, 1) - safeDist];

% Safe Region 3
AE3x = -eye(n+1);
AE3y = [-eye(n+1); eye(n+1)];
AE3z = [-eye(n+1); eye(n+1)];
bE3x = -8*ones(n+1, 1) - safeDist;
bE3y = [zeros(n+1, 1) - safeDist; 9*ones(n+1, 1) - safeDist];
bE3z = [zeros(n+1, 1) - safeDist; 9*ones(n+1, 1) - safeDist];

%% DYNAMIC CONSTRAINTS
Av = [computeD(n, 1); -computeD(n, 1)];
bv = ones(2*n, 1)*velocity_max*segment_time/n;

Aa = [computeD(n, 2); -computeD(n, 2)];
ba = ones(2*(n-1), 1)*acceleration_max*segment_time/(n*(n-1));

%% INITIAL/FINAL CONDITIONS 1
AposInit = [1, zeros(1, n)];
AposFinal = [zeros(1, n), 1];
AvelInit = [-1, 1, zeros(1, n-1)];
AvelFinal = [zeros(1, n-1), -1, 1];
AaccInit = [1, -2, 1, zeros(1, n-2)];
AaccFinal = [zeros(1, n-2), 1, -2, 1];

bposInitx = start_position1(1);
bposInity = start_position1(2);
bposInitz = start_position1(3);
bposFinalx = target_position1(1);
bposFinaly = target_position1(2);
bposFinalz = target_position1(3);
bvelInitx = 0;
bvelInity = 0;
bvelInitz = 0;
bvelFinalx = 0;
bvelFinaly = 0;
bvelFinalz = 0;
baccInitx = 0;
baccInity = 0;
baccInitz = 0;
baccFinalx = 0;
baccFinaly = 0;
baccFinalz = 0;

%% CONTINUITY CONSTRAINTS
AposC = [zeros(1, n), 1, -1, zeros(1, n)];
AvelC = [zeros(1, n-1), -1, 1, 1, -1, zeros(1, n-1)];
AaccC = [zeros(1, n-2), 1, -2, 1, -1, 2, -1, zeros(1, n-2)];

AC = [AposC; AvelC; AaccC];

%% BUILD EQUALITY MATRIX 1
Aeq = [];
for k = 1:segments_num-1
    Aeq = [Aeq, zeros(3*(k-1), n+1); zeros(3, (k-1)*(n+1)), AC];
end
Aeq = [AposInit, zeros(1, (n+1)*(segments_num-1));
       AvelInit, zeros(1, (n+1)*(segments_num-1));
       AaccInit, zeros(1, (n+1)*(segments_num-1));
       Aeq;
       zeros(1, (n+1)*(segments_num-1)), AposFinal;
       zeros(1, (n+1)*(segments_num-1)), AvelFinal;
       zeros(1, (n+1)*(segments_num-1)), AaccFinal;];

Aeq = blkdiag(Aeq, Aeq, Aeq);

beq = [bposInitx; bvelInitx; baccInitx; zeros(3*(segments_num-1), 1); bposFinalx; bvelFinalx; baccFinalx;
       bposInity; bvelInity; baccInity; zeros(3*(segments_num-1), 1); bposFinaly; bvelFinaly; baccFinaly;
       bposInitz; bvelInitz; baccInitz; zeros(3*(segments_num-1), 1); bposFinalz; bvelFinalz; baccFinalz];
% Aeq e beq must be converted to be able to accept d=2 and d=3

%% SAFE REGION CONSTRAINTS
AEx = blkdiag(AE1x, AE2x, AE3x);
AEy = blkdiag(AE1y, AE2y, AE3y);
AEz = blkdiag(AE1z, AE2z, AE3z);

bEx = [bE1x; bE2x; bE3x];
bEy = [bE1y; bE2y; bE3y];
bEz = [bE1z; bE2z; bE3z];

AE = blkdiag(AEx, AEy, AEz);
bE = [bEx; bEy; bEz];

%% BUILD DYNAMIC CONSTRAINT MATRIX

AvFull = [];
for i = 1:segments_num
    AvFull = blkdiag(AvFull, Av);
end

AvFull = [AvFull; -AvFull];
AvFull = blkdiag(AvFull, AvFull, AvFull);

AaFull = [];
for i = 1:segments_num
    AaFull = blkdiag(AaFull, Aa);
end

AaFull = [AaFull; -AaFull];
AaFull = blkdiag(AaFull, AaFull, AaFull);

Adyn = [AvFull; AaFull];

bdyn = [repmat(bv, 2*d*segments_num, 1); repmat(ba, 2*d*segments_num, 1)];

%% FIND INITIAL TRAJECTORY 1
% Cost
H = computeH(n, r);
D = computeD(n, r);
DHD = D'*H*D;
Hcost = [];
for k = 1:segments_num*d
    Hcost = blkdiag(Hcost, DHD);
end

x = quadprog(Hcost, [], [Adyn; AE], [bdyn;bE], Aeq, beq);

%% STORE INITIAL TRAJECTORY 1

trajectory1 = struct;
for i = 1:segments_num
    trajectory1(i).points = zeros(d, n+1);
    for j = 1:d
    trajectory1(i).points(j, :) = x((j-1)*(n+1)*segments_num + (i-1)*(n+1)+1: ...
                                   (j-1)*(n+1)*segments_num + i*(n+1));
    end
    trajectory1(i).time = segment_time;
    trajectory1(i).curve = zeros(d, length(t));
    
    for tt = 1:length(t)
        trajectory1(i).curve(:, tt) = zeros(d, 1);
        for k = 1:n+1
            trajectory1(i).curve(:, tt) = trajectory1(i).curve(:, tt) + bernsteinPol(n, k-1, t(tt))*trajectory1(i).points(:, k);
        end
    end
end

%% PLOTS

figure(1)
plot3(start_position1(1), start_position1(2), start_position1(3), '.g', 'MarkerSize', 50);
hold on
plot3(target_position1(1), target_position1(2), target_position1(3), '.r', 'MarkerSize', 50);
plotWall();

for i = 1:segments_num
    plot3(trajectory1(i).points(1, :), trajectory1(i).points(2, :), trajectory1(i).points(3, :), 'o', 'MarkerSize', 10);
    plot3(trajectory1(i).curve(1, :), trajectory1(i).curve(2, :), trajectory1(i).curve(3, :), 'LineWidth', 5);
%     for tt = 1:10:length(t)
%         plotSphere(trajectory(i).curve(1, tt), trajectory(i).curve(2, tt), trajectory(i).curve(3, tt), safeDist)
%     end
end
hold off

xlim([-5, 20])
ylim([0, 9])
zlim([0, 9])

%% INITIAL/FINAL CONDITIONS 2
AposInit = [1, zeros(1, n)];
AposFinal = [zeros(1, n), 1];
AvelInit = [-1, 1, zeros(1, n-1)];
AvelFinal = [zeros(1, n-1), -1, 1];
AaccInit = [1, -2, 1, zeros(1, n-2)];
AaccFinal = [zeros(1, n-2), 1, -2, 1];

bposInitx = start_position2(1);
bposInity = start_position2(2);
bposInitz = start_position2(3);
bposFinalx = target_position2(1);
bposFinaly = target_position2(2);
bposFinalz = target_position2(3);
bvelInitx = 0;
bvelInity = 0;
bvelInitz = 0;
bvelFinalx = 0;
bvelFinaly = 0;
bvelFinalz = 0;
baccInitx = 0;
baccInity = 0;
baccInitz = 0;
baccFinalx = 0;
baccFinaly = 0;
baccFinalz = 0;

%% BUILD EQUALITY MATRIX 2
Aeq = [];
for k = 1:segments_num-1
    Aeq = [Aeq, zeros(3*(k-1), n+1); zeros(3, (k-1)*(n+1)), AC];
end
Aeq = [AposInit, zeros(1, (n+1)*(segments_num-1));
       AvelInit, zeros(1, (n+1)*(segments_num-1));
       AaccInit, zeros(1, (n+1)*(segments_num-1));
       Aeq;
       zeros(1, (n+1)*(segments_num-1)), AposFinal;
       zeros(1, (n+1)*(segments_num-1)), AvelFinal;
       zeros(1, (n+1)*(segments_num-1)), AaccFinal;];

Aeq = blkdiag(Aeq, Aeq, Aeq);

beq = [bposInitx; bvelInitx; baccInitx; zeros(3*(segments_num-1), 1); bposFinalx; bvelFinalx; baccFinalx;
       bposInity; bvelInity; baccInity; zeros(3*(segments_num-1), 1); bposFinaly; bvelFinaly; baccFinaly;
       bposInitz; bvelInitz; baccInitz; zeros(3*(segments_num-1), 1); bposFinalz; bvelFinalz; baccFinalz];
% Aeq e beq must be converted to be able to accept d=2 and d=3

%% FIND INITIAL TRAJECTORY 2
% Cost
H = computeH(n, r);
D = computeD(n, r);
DHD = D'*H*D;
Hcost = [];
for k = 1:segments_num*d
    Hcost = blkdiag(Hcost, DHD);
end

x = quadprog(Hcost, [], [Adyn; AE], [bdyn;bE], Aeq, beq);

%% STORE INITIAL TRAJECTORY 2

trajectory2 = struct;
for i = 1:segments_num
    trajectory2(i).points = zeros(d, n+1);
    for j = 1:d
    trajectory2(i).points(j, :) = x((j-1)*(n+1)*segments_num + (i-1)*(n+1)+1: ...
                                   (j-1)*(n+1)*segments_num + i*(n+1));
    end
    trajectory2(i).time = segment_time;
    trajectory2(i).curve = zeros(d, length(t));
    
    for tt = 1:length(t)
        trajectory2(i).curve(:, tt) = zeros(d, 1);
        for k = 1:n+1
            trajectory2(i).curve(:, tt) = trajectory2(i).curve(:, tt) + bernsteinPol(n, k-1, t(tt))*trajectory2(i).points(:, k);
        end
    end
end

%% PLOTS

figure(1)
plot3(start_position2(1), start_position2(2), start_position2(3), '.g', 'MarkerSize', 50);
hold on
plot3(target_position2(1), target_position2(2), target_position2(3), '.r', 'MarkerSize', 50);
plotWall();

for i = 1:segments_num
    plot3(trajectory2(i).points(1, :), trajectory2(i).points(2, :), trajectory2(i).points(3, :), 'o', 'MarkerSize', 10);
    plot3(trajectory2(i).curve(1, :), trajectory2(i).curve(2, :), trajectory2(i).curve(3, :), 'LineWidth', 5);
%     for tt = 1:10:length(t)
%         plotSphere(trajectory(i).curve(1, tt), trajectory(i).curve(2, tt), trajectory(i).curve(3, tt), safeDist)
%     end
end
hold off

xlim([-5, 20])
ylim([0, 9])
zlim([0, 9])

%% OBSTACLE PATH
obstaclePoints = [  0,   2,   4,   6,   9,  10;
                    6,   4,   4,   4,   4,   3;
                    9,   7,   5,   5,   5,   3];

obstacleCurve = zeros(d, length(t));
for tt = 1:length(t)
    obstacleCurve(:, tt) = zeros(d, 1);
    for k = 1:n+1
        obstacleCurve(:, tt) = obstacleCurve(:, tt) + bernsteinPol(n, k-1, t(tt))*obstaclePoints(:, k);
    end
end

%% PLOTS

figure(2)
plot3(start_position1(1), start_position1(2), start_position1(2), '.g', 'MarkerSize', 50);
hold on
plot3(target_position1(1), target_position1(2), target_position1(3), '.r', 'MarkerSize', 50);
plot3(start_position2(1), start_position2(2), start_position2(3), '.b', 'MarkerSize', 50);
plot3(target_position2(1), target_position2(2), target_position2(3), '.k', 'MarkerSize', 50);
plotWall();

for i = 1:segments_num
    plot3(trajectory1(i).points(1, :), trajectory1(i).points(2, :), trajectory1(i).points(3, :), 'o', 'MarkerSize', 10);
    plot3(trajectory1(i).curve(1, :), trajectory1(i).curve(2, :), trajectory1(i).curve(3, :), 'LineWidth', 5);
    for tt = 1:10:length(t)
        plotSphere(trajectory1(i).curve(1, tt), trajectory1(i).curve(2, tt), trajectory1(i).curve(3, tt), safeDist)
    end
end
plot3(obstacleCurve(1, :), obstacleCurve(2, :), obstacleCurve(3, :), 'LineWidth', 5);
hold off

xlim([-5, 20])
ylim([0, 9])
zlim([0, 9])

%% %% OPTIMIZATION P/U (ONLY U?)
idx_collision = 2;
oldPoints = trajectory1(idx_collision).points;


Udrone1 = [0, 0.3, 0.7, 1];
Uobstacle = [0, 0.3, 0.6, 1];

Q = computeQ(n*m);

J = @(x) costU(x, m);
Co = BezierComposition(obstaclePoints(:, 1), obstaclePoints(:, 2), ...
                       obstaclePoints(:, 3), obstaclePoints(:, 4), ...
                       obstaclePoints(:, 5), obstaclePoints(:, 6), ...
                       Uobstacle(1), Uobstacle(2), ...
                       Uobstacle(3), Uobstacle(4));
Co = reshape(Co, d, []);
constr = @(x) nlCnstU(x, oldPoints, Co, n, m , d, safeDist, Q);


Aeq = [1, zeros(1, m); zeros(1, m), 1];
beq = [0; 1];

lb = zeros(m+1, 1);
ub = ones(m+1, 1);

% Add dynamic contraints

opts = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluation', 1e4, 'MaxIterations', 1e4);
tic
x = fmincon(J, Udrone1, [], [], Aeq, beq, lb, ub, constr, opts);
toc
nlCnstU(x, oldPoints, Co, n, m , d, safeDist, Q);
UNew = x;

%% COMPUTE DISTANCE
C = BezierComposition(oldPoints(:, 1), oldPoints(:, 2), ...
                      oldPoints(:, 3), oldPoints(:, 4), ...
                      oldPoints(:, 5), oldPoints(:, 6), ...
                      UNew(1), UNew(2), UNew(3), UNew(4));
C = reshape(C, d, []);
diff = C - Co;

normPoints = zeros(2*n*m+1, 1);
for k = 0:2*n*m
    for j = 1:d
        normPoints(k+1) = normPoints(k+1) + diff(j, :)*Q(:, :, k+1)*diff(j, :)';
    end
end

normBez = zeros(1, length(t));
for tt = 1:length(t)
    normBez(tt) = 0;
    for k = 1:2*(n*m)+1
        normBez(tt) = normBez(tt) + bernsteinPol(2*n*m, k-1, t(tt))*normPoints(k);
    end
end

%% PLOT DISTANCE

figure(19)
plot(normBez);
hold on
plot(safeDist^2*ones(length(normBez), 1))
hold off

%% CURVES
bezierCurveD = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveD(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*UNew(k);
    end
    for k = 1:n+1
        bezierCurveD(:, tt) = bezierCurveD(:, tt) + bernsteinPol(n, k-1, u)*oldPoints(:, k);
    end
end

bezierCurveO = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveO(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Uobstacle(k);
    end
    for k = 1:n+1
        bezierCurveO(:, tt) = bezierCurveO(:, tt) + bernsteinPol(n, k-1, u)*obstaclePoints(:, k);
    end
end

% for tt = 1:length(t)
%     figure(20)
%     plot3(bezierCurveD(1, :), bezierCurveD(2, :), bezierCurveD(3, :));
%     hold on
%     plot3(bezierCurveO(1, :), bezierCurveO(2, :), bezierCurveO(3, :));
%     
%     plot3(bezierCurveD(1, tt), bezierCurveD(2, tt), bezierCurveD(3, tt), 'x', 'MarkerSize', 5)
%     plot3(bezierCurveO(1, tt), bezierCurveO(2, tt), bezierCurveO(3, tt), 'x', 'MarkerSize', 5)
%     plotWall()
%     plotSphere(bezierCurveD(1, tt), bezierCurveD(2, tt), bezierCurveD(3, tt), safeDist)
%     hold off
%     xlim([-5, 20])
%     ylim([0, 9])
%     zlim([0, 9])
%     drawnow
% end