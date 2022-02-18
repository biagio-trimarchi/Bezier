%% %% %% AVOIDANCE CLUTTERED %% %% %%
% Notes

%% SETUP
clear 
clc 
addpath('./Functions')

%% PARAMETERS
n = 5;              % Spatial Beziér Curve Order
m = 3;              % Temporal Beziér Curve Order
d = 2;              % Space Dimension
r = 2;              % Shape Parameter
segments_num = 5;   % Number of Segments
t = 0:0.01:1;       % Time Vector (for visualization)

%% %% INITIAL MISSION PLANNING
%% MISSION PARAMETERS
start_position = [1; 4];                    % Start Position
target_position = [8; 4];                   % Target Position
segment_time = 2;                           % Time Allocated to each Segment
total_time = segments_num*segment_time;     % Total Mission Time

%% VEHICLE PARAMETERS
velocity_max = 1000;          % Maximum Velocity
acceleration_max = 100;       % Maximum Acceleration
safeDist = 0.25;            % SafeDistance

%% SAFE REGIONS
% Safe Region 1
AE1x = eye(n+1);
AE1y = -eye(n+1);
bE1x = 2*ones(n+1, 1) - safeDist;
bE1y = zeros(n+1, 1) - safeDist;

% Safe Region 2
AE2x = eye(n+1);
AE2y = [-eye(n+1); eye(n+1)];
bE2x = 5*ones(n+1, 1) - safeDist;
bE2y = [-5*ones(n+1, 1) - safeDist; 6*ones(n+1, 1) - safeDist];

% Safe Region 3
AE3x = [-eye(n+1); eye(n+1)];
AE3y = -eye(n+1);
bE3x = [-4*ones(n+1, 1)-safeDist; 5*ones(n+1, 1)-safeDist];
bE3y = zeros(n+1, 1)-safeDist;

% Safe Region 4
AE4x = -eye(n+1);
AE4y = [-eye(n+1); eye(n+1)];
bE4x = -4*ones(n+1, 1)-safeDist;
bE4y = [zeros(n+1, 1)-safeDist; ones(n+1, 1)-safeDist];

% Safe Region 5
AE5x = -eye(n+1);
AE5y = -eye(n+1);
bE5x = -7*ones(n+1, 1)-safeDist;
bE5y = zeros(n+1, 1)-safeDist;

%% DYNAMIC CONSTRAINTS
Av = [computeD(n, 1); -computeD(n, 1)];
bv = ones(2*n, 1)*velocity_max*segment_time/n;

Aa = [computeD(n, 2); -computeD(n, 2)];
ba = ones(2*(n-1), 1)*acceleration_max*segment_time/(n*(n-1));

%% INITIAL/FINAL CONDITIONS
AposInit = [1, zeros(1, n)];
AposFinal = [zeros(1, n), 1];
AvelInit = [-1, 1, zeros(1, n-1)];
AvelFinal = [zeros(1, n-1), -1, 1];
AaccInit = [1, -2, 1, zeros(1, n-2)];
AaccFinal = [zeros(1, n-2), 1, -2, 1];

bposInitx = start_position(1);
bposInity = start_position(2);
bposFinalx = target_position(1);
bposFinaly = target_position(2);
bvelInitx = 0;
bvelInity = 0;
bvelFinalx = 0;
bvelFinaly = 0;
baccInitx = 0;
baccInity = 0;
baccFinalx = 0;
baccFinaly = 0;

%% CONTINUITY CONSTRAINTS
AposC = [zeros(1, n), 1, -1, zeros(1, n)];
AvelC = [zeros(1, n-1), -1, 1, 1, -1, zeros(1, n-1)];
AaccC = [zeros(1, n-2), 1, -2, 1, -1, 2, -1, zeros(1, n-2)];

AC = [AposC; AvelC; AaccC];

%% BUILD EQUALITY MATRIX
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

Aeq = blkdiag(Aeq, Aeq);

beq = [bposInitx; bvelInitx; baccInitx; zeros(3*(segments_num-1), 1); bposFinalx; bvelFinalx; baccFinalx;
       bposInity; bvelInity; baccInity; zeros(3*(segments_num-1), 1); bposFinaly; bvelFinaly; baccFinaly];
% Aeq e beq must be converted to be able to accept d=2 and d=3
   
%% SAFE REGION CONSTRAINTS
AEx = blkdiag(AE1x, AE2x, AE3x, AE4x, AE5x);
AEy = blkdiag(AE1y, AE2y, AE3y, AE4y, AE5y);

bEx = [bE1x; bE2x; bE3x; bE4x; bE5x];
bEy = [bE1y; bE2y; bE3y; bE4y; bE5y];

AE = blkdiag(AEx, AEy);
bE = [bEx; bEy];

%% BUILD DYNAMIC CONSTRAINT MATRIX

AvFull = [];
for i = 1:segments_num
    AvFull = blkdiag(AvFull, Av);
end

AvFull = [AvFull; -AvFull];
AvFull = blkdiag(AvFull, AvFull);

AaFull = [];
for i = 1:segments_num
    AaFull = blkdiag(AaFull, Aa);
end

AaFull = [AaFull; -AaFull];
AaFull = blkdiag(AaFull, AaFull);

Adyn = [AvFull; AaFull];

bdyn = [repmat(bv, 2*d*segments_num, 1); repmat(ba, 2*d*segments_num, 1)];

%% FIND INITIAL TRAJECTORY
% Cost
H = computeH(n, r);
D = computeD(n, r);
DHD = D'*H*D;
Hcost = [];
for k = 1:segments_num*d
    Hcost = blkdiag(Hcost, DHD);
end

x = quadprog(Hcost, [], [Adyn; AE], [bdyn;bE], Aeq, beq);

%% STORE INITIAL TRAJECTORY

trajectory = struct;
for i = 1:segments_num
    trajectory(i).points = zeros(d, n+1);
    for j = 1:d
    trajectory(i).points(j, :) = x((j-1)*(n+1)*segments_num + (i-1)*(n+1)+1: ...
                                   (j-1)*(n+1)*segments_num + i*(n+1));
    end
    trajectory(i).time = segment_time;
    trajectory(i).curve = zeros(d, length(t));
    
    for tt = 1:length(t)
        trajectory(i).curve(:, tt) = zeros(d, 1);
        for k = 1:n+1
            trajectory(i).curve(:, tt) = trajectory(i).curve(:, tt) + bernsteinPol(n, k-1, t(tt))*trajectory(i).points(:, k);
        end
    end
end


%% PLOTS

figure(1)
plot(start_position(1), start_position(2), 'x', 'MarkerSize', 10);
hold on
plot(target_position(1), target_position(2), 'x', 'MarkerSize', 10);
drawObstacles();

for i = 1:segments_num
    plot(trajectory(i).points(1, :), trajectory(i).points(2, :), 'o', 'MarkerSize', 10);
    plot(trajectory(i).curve(1, :), trajectory(i).curve(2, :));
end
hold off

xlim([0, 9])
ylim([0, 9])

%% %% OBSTACLE AVOIDANCE CLASSIC

%% OBSTACLE PATH
obstaclePoints = [4.5, 4.5, 4.5, 4.5, 4.5, 4.5;
                    8,   3,   3,   3,   2,   2];

obstacleCurve = zeros(d, length(t));
for tt = 1:length(t)
    obstacleCurve(:, tt) = zeros(d, 1);
    for k = 1:n+1
        obstacleCurve(:, tt) = obstacleCurve(:, tt) + bernsteinPol(n, k-1, t(tt))*obstaclePoints(:, k);
    end
end


%% PLOT

figure(2)
plot(start_position(1), start_position(2), 'x', 'MarkerSize', 10);
hold on
plot(target_position(1), target_position(2), 'x', 'MarkerSize', 10);
drawObstacles();

for i = 1:segments_num
    plot(trajectory(i).points(1, :), trajectory(i).points(2, :), 'o', 'MarkerSize', 10);
    plot(trajectory(i).curve(1, :), trajectory(i).curve(2, :));
end
plot(obstacleCurve(1, :), obstacleCurve(2, :));
hold off

xlim([0, 9])
ylim([0, 9])

%% PLOT ANIMATION

for tt = 1:length(t)
    figure(3)
    plot(trajectory(3).curve(1, :), trajectory(3).curve(2, :));
    hold on
    plot(obstacleCurve(1, :), obstacleCurve(2, :));
    
    plot(trajectory(3).curve(1, tt), trajectory(3).curve(2, tt), 'o', 'MarkerSize', 10)
    plot(obstacleCurve(1, tt), obstacleCurve(2, tt), 'o', 'MarkerSize', 10)
    drawObstacles()
    plotCircle(trajectory(3).curve(1, tt), trajectory(3).curve(2, tt), safeDist)
    hold off
    xlim([0, 9])
    ylim([0, 9])
    drawnow
end


%% DIVIDE CURVE AGENT

tau = 0.8;          % Moment of collision
idx_collision = 3;

[Px1, Px2] = bezierDivision(trajectory(idx_collision).points(1, :), tau, n);
[Py1, Py2] = bezierDivision(trajectory(idx_collision).points(2, :), tau, n);

for i=segments_num:-1:idx_collision+1
    trajectory(i+1) = trajectory(i);
end
segments_num = segments_num + 1;

trajectory(idx_collision).points(1, :) = Px1;
trajectory(idx_collision).points(2, :) = Py1;
trajectory(idx_collision).time = segment_time*tau;
trajectory(idx_collision).curve = zeros(d, length(t));

for tt = 1:length(t)
    trajectory(idx_collision).curve(:, tt) = zeros(d, 1);
    for k = 1:n+1
        trajectory(idx_collision).curve(:, tt) = trajectory(idx_collision).curve(:, tt) + bernsteinPol(n, k-1, t(tt))*trajectory(idx_collision).points(:, k);
    end
end

trajectory(idx_collision+1).points(1, :) = Px2;
trajectory(idx_collision+1).points(2, :) = Py2;
trajectory(idx_collision+1).time = segment_time*(1-tau);
trajectory(idx_collision+1).curve = zeros(d, length(t));

for tt = 1:length(t)
    trajectory(idx_collision+1).curve(:, tt) = zeros(d, 1);
    for k = 1:n+1
        trajectory(idx_collision+1).curve(:, tt) = trajectory(idx_collision+1).curve(:, tt) + bernsteinPol(n, k-1, t(tt))*trajectory(idx_collision+1).points(:, k);
    end
end

%% DIVIDE CURVE OBSTACLE

[Pox1, Pox2] = bezierDivision(obstaclePoints(1, :), tau, n);
[Poy1, Poy2] = bezierDivision(obstaclePoints(2, :), tau, n);


obstaclePointsDiv1 = [Pox1'; Poy1'];
obstaclePointsDiv2 = [Pox2'; Poy2'];

obstacleCurveDiv1 = zeros(d, length(t));
for tt = 1:length(t)
    obstacleCurveDiv1(:, tt) = zeros(d, 1);
    for k = 1:n+1
        obstacleCurveDiv1(:, tt) = obstacleCurveDiv1(:, tt) + bernsteinPol(n, k-1, t(tt))*obstaclePointsDiv1(:, k);
    end
end

obstacleCurveDiv2 = zeros(d, length(t));
for tt = 1:length(t)
    obstacleCurveDiv2(:, tt) = zeros(d, 1);
    for k = 1:n+1
        obstacleCurveDiv2(:, tt) = obstacleCurveDiv2(:, tt) + bernsteinPol(n, k-1, t(tt))*obstaclePointsDiv2(:, k);
    end
end

%% PLOT

figure(4)
plot(start_position(1), start_position(2), 'x', 'MarkerSize', 10);
hold on
plot(target_position(1), target_position(2), 'x', 'MarkerSize', 10);
drawObstacles();

for i = 1:segments_num
    plot(trajectory(i).points(1, :), trajectory(i).points(2, :), 'o', 'MarkerSize', 10);
    plot(trajectory(i).curve(1, :), trajectory(i).curve(2, :));
end
plot(obstacleCurveDiv1(1, :), obstacleCurveDiv1(2, :));
plot(obstacleCurveDiv2(1, :), obstacleCurveDiv2(2, :));
hold off

xlim([0, 9])
ylim([0, 9])

%% SETUP OPTIMIZATION PROBLEM

segments__opt_num = 2;
Hcost = [];
for k = 1:segments__opt_num*d
    Hcost = blkdiag(Hcost, DHD);
end

J = @(x) x'*Hcost*x;
%% EQUALITY CONTRAINT MATRIX
Aeq = [];
for k = 1:segments__opt_num-1
    Aeq = [Aeq, zeros(3*(k-1), n+1); zeros(3, (k-1)*(n+1)), AC];
end
Aeq = [AposInit, zeros(1, (n+1)*(segments__opt_num-1));
       AvelInit, zeros(1, (n+1)*(segments__opt_num-1));
       AaccInit, zeros(1, (n+1)*(segments__opt_num-1));
       Aeq;
       zeros(1, (n+1)*(segments__opt_num-1)), AposFinal;
       zeros(1, (n+1)*(segments__opt_num-1)), AvelFinal;
       zeros(1, (n+1)*(segments__opt_num-1)), AaccFinal;];

Aeq = blkdiag(Aeq, Aeq);

posInit = trajectory(idx_collision).points(:, 1);
velInit = trajectory(idx_collision).points(:, 2) - trajectory(idx_collision).points(:, 1);
accInit = trajectory(idx_collision).points(:, 3) - 2*trajectory(idx_collision).points(:, 2) + trajectory(idx_collision).points(:, 1);

posFinal = trajectory(idx_collision+1).points(:, end);
velFinal = trajectory(idx_collision+1).points(:, end) - trajectory(idx_collision).points(:, end-1);
accFinal = trajectory(idx_collision+1).points(:, end) - 2*trajectory(idx_collision).points(:, end-1) + trajectory(idx_collision).points(:, end-2);

beq = [posInit(1); velInit(1); accInit(1); zeros(3*(segments__opt_num-1), 1); posFinal(1); velFinal(1); accFinal(1);
       posInit(2); velInit(2); accInit(2); zeros(3*(segments__opt_num-1), 1); posFinal(2); velFinal(2); accFinal(2)];
% Aeq e beq must be converted to be able to accept d=2 and d=3

%% SAFE REGION CONSTRAINTS
AEx = blkdiag(AE3x, AE3x);
AEy = blkdiag(AE3y, AE3y);

bEx = [bE3x; bE3x];
bEy = [bE3y; bE3y];

AE = blkdiag(AEx, AEy);
bE = [bEx; bEy];

%% DYNAMIC CONSTRAINTS
AvFull = [];
for i = 1:segments__opt_num
    AvFull = blkdiag(AvFull, Av);
end

AvFull = [AvFull; -AvFull];
AvFull = blkdiag(AvFull, AvFull);

AaFull = [];
for i = 1:segments__opt_num
    AaFull = blkdiag(AaFull, Aa);
end

AaFull = [AaFull; -AaFull];
AaFull = blkdiag(AaFull, AaFull);

Adyn = [AvFull; AaFull];

bdyn = [repmat(bv, 2*d*segments__opt_num, 1); repmat(ba, 2*d*segments__opt_num, 1)];

%% REPLAN
x0 = [trajectory(idx_collision).points(1, :), trajectory(idx_collision+1).points(1, :), trajectory(idx_collision).points(2, :), trajectory(idx_collision+1).points(2, :)]';
opts = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluation', inf, 'MaxIterations', inf);

Q = computeQ(n);
constr = @(x) nlCnstP(x, obstaclePointsDiv1, obstaclePointsDiv2, n, d, safeDist, Q);

tic
[x, fval] = fmincon(J, x0, [Adyn; AE], [bdyn;bE], Aeq, beq, [], [], constr, opts);
toc
x1 = quadprog(Hcost, [], [Adyn; AE], [bdyn;bE], Aeq, beq);


%% STORE REPLANNED PATH
for i = 1:segments__opt_num
    for j = 1:d
    trajectory(i+idx_collision-1).points(j, :) = x((j-1)*(n+1)*segments__opt_num + (i-1)*(n+1)+1: ...
                                   (j-1)*(n+1)*segments__opt_num + i*(n+1));
    end

    trajectory(i+idx_collision-1).curve = zeros(d, length(t)); 
    for tt = 1:length(t)
        trajectory(i+idx_collision-1).curve(:, tt) = zeros(d, 1);
        for k = 1:n+1
            trajectory(i+idx_collision-1).curve(:, tt) = trajectory(i+idx_collision-1).curve(:, tt) + bernsteinPol(n, k-1, t(tt))*trajectory(i+idx_collision-1).points(:, k);
        end
    end
end

%% PLOTS
for tt = 1:length(t)
    figure(10)
    plot(trajectory(3).curve(1, :), trajectory(3).curve(2, :));
    hold on
    plot(obstacleCurveDiv1(1, :), obstacleCurveDiv1(2, :));
    
    plot(trajectory(3).curve(1, tt), trajectory(3).curve(2, tt), 'o', 'MarkerSize', 10)
    plot(obstacleCurveDiv1(1, tt), obstacleCurveDiv1(2, tt), 'o', 'MarkerSize', 10)
    drawObstacles()
    plotCircle(trajectory(3).curve(1, tt), trajectory(3).curve(2, tt), safeDist)
    hold off
    xlim([0, 9])
    ylim([0, 9])
    drawnow
end

for tt = 1:length(t)
    figure(11)
    plot(trajectory(4).curve(1, :), trajectory(4).curve(2, :));
    hold on
    plot(obstacleCurveDiv2(1, :), obstacleCurveDiv2(2, :));
    
    plot(trajectory(4).curve(1, tt), trajectory(4).curve(2, tt), 'o', 'MarkerSize', 10)
    plot(obstacleCurveDiv2(1, tt), obstacleCurveDiv2(2, tt), 'o', 'MarkerSize', 10)
    drawObstacles()
    plotCircle(trajectory(4).curve(1, tt), trajectory(4).curve(2, tt), safeDist)
    hold off
    xlim([0, 9])
    ylim([0, 9])
    drawnow
end

figure(20)
plot(start_position(1), start_position(2), 'x', 'MarkerSize', 10);
hold on
plot(target_position(1), target_position(2), 'x', 'MarkerSize', 10);
drawObstacles();

for i = 1:segments_num
    plot(trajectory(i).points(1, :), trajectory(i).points(2, :), 'o', 'MarkerSize', 10);
    plot(trajectory(i).curve(1, :), trajectory(i).curve(2, :));
end
hold off

xlim([0, 9])
ylim([0, 9])

%% %% OPTIMIZATION P/U (ONLY U?)

