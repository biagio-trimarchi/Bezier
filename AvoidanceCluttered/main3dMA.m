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
start_position2 = [0; 2; 7];                 % Start Position
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

figure(2)
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


%% PLOTS

figure(3)
plot3(start_position1(1), start_position1(2), start_position1(2), '.g', 'MarkerSize', 50);
hold on
plot3(start_position2(1), start_position2(2), start_position2(3), '.g', 'MarkerSize', 50);
plot3(target_position1(1), target_position1(2), target_position1(3), '.r', 'MarkerSize', 50);
plot3(target_position2(1), target_position2(2), target_position2(3), '.r', 'MarkerSize', 50);
plotWall();

for i = 1:segments_num
    plot3(trajectory1(i).points(1, :), trajectory1(i).points(2, :), trajectory1(i).points(3, :), 'o', 'MarkerSize', 10);
    plot3(trajectory1(i).curve(1, :), trajectory1(i).curve(2, :), trajectory1(i).curve(3, :), 'LineWidth', 5);
    plot3(trajectory2(i).points(1, :), trajectory2(i).points(2, :), trajectory2(i).points(3, :), 'o', 'MarkerSize', 10);
    plot3(trajectory2(i).curve(1, :), trajectory2(i).curve(2, :), trajectory2(i).curve(3, :), 'LineWidth', 5);
    for tt = 1:10:length(t)
        plotSphere(trajectory1(i).curve(1, tt), trajectory1(i).curve(2, tt), trajectory1(i).curve(3, tt), safeDist)
    end
end
hold off

xlim([-5, 20])
ylim([0, 9])
zlim([0, 9])

%% %% OPTIMIZATION P/U (ONLY U?)
idx_collision = 2;
oldPoints1 = trajectory1(idx_collision).points;
oldPoints2 = trajectory2(idx_collision).points;

Udrone1 = [0, 0.3, 0.7, 1];
Udrone2 = [0, 0.2, 0.4, 1];

Q = computeQ(n*m);

J = @(x) costU(x(1:m+1), m) + costU(x(m+2:end), m);
constr = @(x) nlCnstU2(x, oldPoints1, oldPoints2, n, m , d, safeDist, Q);


Aeq = [1, zeros(1, m); zeros(1, m), 1];
Aeq = blkdiag(Aeq, Aeq);
beq = [0; 1; 0; 1];

lb = zeros(2*(m+1), 1);
ub = ones(2*(m+1), 1);

% Add dynamic contraints

opts = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluation', 1e4, 'MaxIterations', 1e4);
tic
x = fmincon(J, [Udrone1, Udrone2], [], [], Aeq, beq, lb, ub, constr, opts);
toc
UNew1 = x(1:m+1);
UNew2 = x(m+2:end);

%% COMPUTE DISTANCE BEFORE
C1 = BezierComposition(oldPoints1(:, 1), oldPoints1(:, 2), ...
                      oldPoints1(:, 3), oldPoints1(:, 4), ...
                      oldPoints1(:, 5), oldPoints1(:, 6), ...
                      Udrone1(1), Udrone1(2), Udrone1(3), Udrone1(4));
C1 = reshape(C1, d, []);

C2 = BezierComposition(oldPoints2(:, 1), oldPoints2(:, 2), ...
                      oldPoints2(:, 3), oldPoints2(:, 4), ...
                      oldPoints2(:, 5), oldPoints2(:, 6), ...
                      Udrone2(1), Udrone2(2), Udrone2(3), Udrone2(4));
C2 = reshape(C2, d, []);
diff = C1 - C2;

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

%% PLOT DISTANCE BEFORE

figure(18)
plot(normBez);
hold on
plot(safeDist^2*ones(length(normBez), 1))
hold off

%% COMPUTE DISTANCE AFTER
C1 = BezierComposition(oldPoints1(:, 1), oldPoints1(:, 2), ...
                      oldPoints1(:, 3), oldPoints1(:, 4), ...
                      oldPoints1(:, 5), oldPoints1(:, 6), ...
                      UNew1(1), UNew1(2), UNew1(3), UNew1(4));
C1 = reshape(C1, d, []);

C2 = BezierComposition(oldPoints2(:, 1), oldPoints2(:, 2), ...
                      oldPoints2(:, 3), oldPoints2(:, 4), ...
                      oldPoints2(:, 5), oldPoints2(:, 6), ...
                      UNew2(1), UNew2(2), UNew2(3), UNew2(4));
C2 = reshape(C2, d, []);
diff = C1 - C2;

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

%% PLOT DISTANCE AFTER

figure(19)
plot(normBez);
hold on
plot(safeDist^2*ones(length(normBez), 1))
hold off

%% CURVES
bezierCurveD1 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveD1(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*UNew1(k);
    end
    for k = 1:n+1
        bezierCurveD1(:, tt) = bezierCurveD1(:, tt) + bernsteinPol(n, k-1, u)*oldPoints1(:, k);
    end
end

bezierCurveD2 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurveD2(:, tt) = zeros(d, 1);
    u = 0;
    for k = 1:m+1
        u = u + bernsteinPol(m, k-1, t(tt))*Udrone2(k);
    end
    for k = 1:n+1
        bezierCurveD2(:, tt) = bezierCurveD2(:, tt) + bernsteinPol(n, k-1, u)*oldPoints2(:, k);
    end
end

% for tt = 1:length(t)
%     figure(20)
%     plot3(bezierCurveD1(1, :), bezierCurveD1(2, :), bezierCurveD1(3, :), 'LineWidth', 5);
%     hold on
%     plot3(bezierCurveD2(1, :), bezierCurveD2(2, :), bezierCurveD2(3, :), 'LineWidth', 5);
%     
%     plot3(bezierCurveD1(1, tt), bezierCurveD1(2, tt), bezierCurveD1(3, tt), 'x', 'MarkerSize', 20)
%     plot3(bezierCurveD2(1, tt), bezierCurveD2(2, tt), bezierCurveD2(3, tt), 'x', 'MarkerSize', 20)
%     plotWall()
%     plotSphere(bezierCurveD1(1, tt), bezierCurveD1(2, tt), bezierCurveD1(3, tt), safeDist)
%     hold off
%     xlim([-5, 20])
%     ylim([0, 9])
%     zlim([0, 9])
%     drawnow
% end