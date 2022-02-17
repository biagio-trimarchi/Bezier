%% AVOIDANCE CLUTTERED
% Notes

%% SETUP
clear 
clc 
addpath('./Functions')

%% PARAMETERS
n = 5;              % Spatial Beziér Curve Order
m = 3;              % Temporal Beziér Curve Order
d = 2;              % Space Dimension
r = 5;              % Shape Parameter
segments_num = 5;   % Number of Segments
t = 0:0.01:1;       % Time Vector (for visualization)

%% %% INITIAL MISSION PLANNING
%% MISSION PARAMETERS
start_position = [1; 4];                    % Start Position
target_position = [8; 4];                   % Target Position
segment_time = 2;                           % Time Allocated to each Segment
total_time = segments_num*segment_time;     % Total Mission Time

%% VEHICLE PARAMETERS
velocity_max = 10;          % Maximum Velocity
acceleration_max = 5;       % Maximum Acceleration
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

%% Dynamic Constraints
Av = [computeD(n, 1); -computeD(n, 1)];
bv = ones(2*n, 1)*velocity_max*segment_time/n;

Aa = [computeD(n, 2); -computeD(n, 2)];
ba = ones(2*(n-1), 1)*acceleration_max*segment_time/(n*(n-1));

%% Initial/Final Position
AposInit = [1, zeros(1, n)];
AposFinal = [zeros(1, n), 1];

bposInitx = [start_position(1)];
bposInity = [start_position(2)];
bposFinalx = [target_position(1)];
bposFinaly = [target_position(2)];

%% Continuity Constraints
AposC = [zeros(1, n), 1, -1, zeros(1, n)];
AvelC = [zeros(1, n-1), -1, 1, 1, -1, zeros(1, n-1)];
AaccC = [zeros(1, n-2), 1, -2, 1, -1, 2, -1, zeros(1, n-2)];

AC = [AposC; AvelC; AaccC];

%% Build Equality Matrix
Aeq = [];
for k = 1:segments_num-1
    Aeq = [Aeq, zeros(3*(k-1), n+1); zeros(3, (k-1)*(n+1)), AC];
end
Aeq = [AposInit, zeros(1, (n+1)*(segments_num-1));
       Aeq;
       zeros(1, (n+1)*(segments_num-1)), AposFinal];

Aeq = blkdiag(Aeq, Aeq);

beq = [bposInitx; zeros(3*(segments_num-1), 1); bposFinalx;
       bposInity; zeros(3*(segments_num-1), 1); bposFinaly];
% Aeq e beq must be converted to be able to accept d=2 and d=3
   
%% SAFE REGION CONSTRAINTS
AEx = blkdiag(AE1x, AE2x, AE3x, AE4x, AE5x);
AEy = blkdiag(AE1y, AE2y, AE3y, AE4y, AE5y);

bEx = [bE1x; bE2x; bE3x; bE4x; bE5x];
bEy = [bE1y; bE2y; bE3y; bE4y; bE5y];

AE = blkdiag(AEx, AEy);
bE = [bEx; bEy];

%% Dynamic Constraints Full

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

%% Find Initial Trajectory
% Cost
H = computeH(n, r);
D = computeD(n, r);
DHD = D'*H*D;
Hcost = [];
for k = 1:segments_num*d
    Hcost = blkdiag(Hcost, DHD);
end

x = quadprog(Hcost, [], [Adyn; AE], [bdyn;bE], Aeq, beq);

%% Store Initial Trajectory

trajectory = struct;
for i = 1:segments_num
    trajectory(i).points = zeros(d, n+1);
    for j = 1:d
    trajectory(i).points(j, :) = [x((j-1)*(n+1)*segments_num + (i-1)*(n+1)+1: ...
                                    (j-1)*(n+1)*segments_num + i*(n+1))];
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
                    8,   7,   3,   2,   1,   1];

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

% for tt = 1:length(t)
%     figure(3)
%     plot(trajectory(3).curve(1, :), trajectory(3).curve(2, :));
%     hold on
%     plot(obstacleCurve(1, :), obstacleCurve(2, :));
%     
%     plot(trajectory(3).curve(1, tt), trajectory(3).curve(2, tt), 'o', 'MarkerSize', 10)
%     plot(obstacleCurve(1, tt), obstacleCurve(2, tt), 'o', 'MarkerSize', 10)
%     drawObstacles()
%     plotCircle(trajectory(3).curve(1, tt), trajectory(3).curve(2, tt), safeDist)
%     hold off
%     xlim([0, 9])
%     ylim([0, 9])
%     drawnow
% end


%% DIVIDE CURVE

tau = 0.5;          % Moment of collision
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
plot(obstacleCurve(1, :), obstacleCurve(2, :));
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
%% Equality Constraints
Aeq = [];
for k = 1:segments__opt_num-1
    Aeq = [Aeq, zeros(3*(k-1), n+1); zeros(3, (k-1)*(n+1)), AC];
end
Aeq = [AposInit, zeros(1, (n+1)*(segments__opt_num-1));
       Aeq;
       zeros(1, (n+1)*(segments__opt_num-1)), AposFinal];

Aeq = blkdiag(Aeq, Aeq);

beq = [trajectory(idx_collision).points(1, 1); zeros(3*(segments__opt_num-1), 1); trajectory(idx_collision+1).points(1, end);
       trajectory(idx_collision).points(2, 1); zeros(3*(segments__opt_num-1), 1); trajectory(idx_collision+1).points(2, end)];
% Aeq e beq must be converted to be able to accept d=2 and d=3

%% Safe Region Constraints
AEx = blkdiag(AE3x, AE3x);
AEy = blkdiag(AE3y, AE3y);

bEx = [bE3x; bE3x];
bEy = [bE3y; bE3y];

AE = blkdiag(AEx, AEy);
bE = [bEx; bEy];

%% Dynamic Contraints
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

%% Replan
x0 = [trajectory(idx_collision).points(1, :), trajectory(idx_collision+1).points(1, :), trajectory(idx_collision).points(2, :), trajectory(idx_collision+1).points(2, :)]';
opts = optimoptions('fmincon','Algorithm','sqp');

[x, fval] = fmincon(J, x0, [Adyn; AE], [bdyn;bE], Aeq, beq);