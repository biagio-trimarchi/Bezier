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
start_position = [-5; 7; 7];                 % Start Position
target_position = [15; 5; 3];               % Target Position
segment_time = 5;                           % Time Allocated to each Segment
total_time = segments_num*segment_time;     % Total Mission Time

%% VEHICLE PARAMETERS
velocity_max = 50;         % Maximum Velocity
acceleration_max = 20;     % Maximum Acceleration
safeDist = 1;              % SafeDistance

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
bE2x = [5*ones(n+1, 1) - safeDist; 8*ones(n+1, 1) - safeDist];
bE2y = [-3*ones(n+1, 1) - safeDist; 6*ones(n+1, 1) - safeDist];
bE2z = [-3*ones(n+1, 1) - safeDist; 6*ones(n+1, 1) - safeDist];

% Safe Region 3
AE3x = -eye(n+1);
AE3y = [-eye(n+1); eye(n+1)];
AE3z = [-eye(n+1); eye(n+1)];
bE3x = -5*ones(n+1, 1) - safeDist;
bE3y = [zeros(n+1, 1) - safeDist; 9*ones(n+1, 1) - safeDist];
bE3z = [zeros(n+1, 1) - safeDist; 9*ones(n+1, 1) - safeDist];

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
bposInitz = start_position(3);
bposFinalx = target_position(1);
bposFinaly = target_position(2);
bposFinalz = target_position(3);
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
plot3(start_position(1), start_position(2), start_position(2), '.g', 'MarkerSize', 50);
hold on
plot3(target_position(1), target_position(2), target_position(3), '.r', 'MarkerSize', 50);
plotWall();

for i = 1:segments_num
    plot3(trajectory(i).points(1, :), trajectory(i).points(2, :), trajectory(i).points(3, :), 'o', 'MarkerSize', 10);
    plot3(trajectory(i).curve(1, :), trajectory(i).curve(2, :), trajectory(i).curve(3, :));
    for tt = 1:10:length(t)
        plotSphere(trajectory(i).curve(1, tt), trajectory(i).curve(2, tt), trajectory(i).curve(3, tt), safeDist)
    end
end
hold off

xlim([-5, 20])
ylim([0, 9])
zlim([0, 9])

