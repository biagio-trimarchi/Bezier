%% INTERPOLATE TRAJECTORIES
clear all;
close all;
clc;

%% Parameters
n = 5;
m = 3;
bag = rosbag('objBagFiles/thirdTest.bag');
selection = select(bag,'Topic','/tf');
msg = readMessages(selection,'DataFormat','struct');
freq = 20;
dt = 1/freq;

%% Read bag file
idx = 1;
for ii = 1:size(msg,1)
    if size(struct2table(msg{ii}.Transforms), 1) == 0
        continue;
    end

    objPose(1, idx) = msg{ii}.Transforms.Transform.Translation.X;
    objPose(2, idx) = msg{ii}.Transforms.Transform.Translation.Y;
    objPose(3, idx) = msg{ii}.Transforms.Transform.Translation.Z;
    objPose(4, idx) = (idx-1)*dt;
    idx = idx + 1;
end

tF = objPose(end, end);
L = 30;                                                    % Number of sampled points
objPose(4, :) = objPose(4,:) / objPose(4, end);            % Scale time between 0 and 1
collectedPoints = objPose(:, 1:L);                         % Select first L points as collected points

%% Least Square Fitting
controlPoints = interpolateBezierLS(collectedPoints, n);    % Compute control points trough least square

%% Compute Bezier Curve
t = 0:0.01:1;
bezierCurve = zeros(m, length(t));
for tt = 1:length(t)
    bezierCurve(:, tt) = zeros(m, 1);
    for k = 1:n+1
        bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(n, k-1, t(tt))*controlPoints(:, k);
    end
end

%% Plots
% Obstacle trajectory
figure(1);
plot3(objPose(1,:), objPose(2,:), objPose(3,:), 'LineWidth', 1.5);
legend('Obstacle')

% Computed Bezier Curve
figure(2)
plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))
legend('Bezier Curve Least Square')

% Bezier Curve and Obstacle
figure(3)
plot3(objPose(1,:), objPose(2,:), objPose(3,:), 'LineWidth', 1.5);
hold on
plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))
plot3(collectedPoints(1, :), collectedPoints(2, :), collectedPoints(3, :), 'o');
hold off
legend('Obstacle', 'Bezier Curve Least Square', 'Collected Points')


%% Fitting with Regularization
tL = collectedPoints(end, end);             % Final instant in sampled points
w = ones(L, 1);                             % Weigths
wp = 0;

% Matrices for regularization cost
[H, D] = buildHD(n, 2, m);                  
J2 = D'*H*D;                                

% Matrices for trajectory fitting
B = buildB(n, collectedPoints(4, :), m, L, w);
Bp = buildBp(n, collectedPoints(1:3, :), collectedPoints(4, :), m, L, w);

% Constraints
A1 = zeros(n-2+1, n+1);
for i = 1:n-2+1
    A1(i,i) = 1;
    A1(i, i+1) = -2;
    A1(i, (i+2)) = 1;
end
A2 = zeros(n-2+1, n+1);
for i = 1:n-2+1
    A2(i,i) = -1;
    A2(i, i+1) = 2;
    A2(i, (i+2)) = -1;
end
A1 = kron(A1, eye(m));
A2 = kron(A2, eye(m));
A = [A1; A2];

b1 = [0.001; 0.001; 0.0];
b1 = kron(ones(n-2+1, 1), b1);
b2 = [0.001; 0.001; 3.0];
b2 = kron(ones(n-2+1, 1), b2);
b = [b1; b2];


% Compute control Points
controlPoints2 = quadprog(B + wp*J2, Bp, A, b);
controlPoints2 = reshape(controlPoints2, 3, []);

%% Bezier Curve
bezierCurve2 = zeros(3, length(t));
for tt = 1:length(t)
    bezierCurve2(:, tt) = [0; 0; 0];
    for k = 1:n+1
        bezierCurve2(:, tt) = bezierCurve2(:, tt) + bernsteinPol(n, k-1, t(tt))*controlPoints2(:, k);
    end
end

%% Plots
% Obstacle trajectory
figure(4);
plot3(objPose(1,:), objPose(2,:), objPose(3,:), 'LineWidth', 1.5);
legend('Obstacle')

% Computed Bezier Curve
figure(5)
plot3(bezierCurve2(1, :), bezierCurve2(2, :), bezierCurve2(3, :))
legend('Bezier Curve Regularization')

% Bezier Curve and Obstacle
figure(6)
plot3(objPose(1,:), objPose(2,:), objPose(3,:), 'LineWidth', 1.5);
hold on
plot3(bezierCurve2(1, :), bezierCurve2(2, :), bezierCurve2(3, :))
plot3(collectedPoints(1, :), collectedPoints(2, :), collectedPoints(3, :), 'o');
hold off
legend('Obstacle', 'Bezier Curve Least Square', 'Collected Points')

%% FUNCTIONS
function controlPoints = interpolateBezierLS(points, n)
    tau = points(end, :);
    points = points(1:end-1, :);
    dim = size(points, 1);
    
    points = reshape(points, [], 1);
    
    B = zeros(size(tau, 2)*dim, (n+1)*dim);
    
    for i=1:size(B,1)/dim
        for j=1:size(B,2)/dim
            auxB = eye(dim)*bernsteinPol(n, j-1, tau(i));
            B(1+(dim*(i-1)):dim*i, 1+(dim*(j-1)):dim*j) = auxB;

        end
    end
    
    % controlPoints = lsqr(B, points);
    controlPoints = pinv(B)*points;
    controlPoints = reshape(controlPoints, dim, []);
end

function controlPoints = interpolateBezierWLS(points, n, W)
    tau = points(end, :);
    points = points(1:end-1, :);
    dim = size(points, 1);
    
    points = reshape(points, [], 1);
    
    B = zeros(size(tau, 2)*dim, (n+1)*dim);
    
    for i=1:size(B,1)/dim
        for j=1:size(B,2)/dim
            auxB = eye(dim)*bernsteinPol(n, j-1, tau(i));
            B(1+(dim*(i-1)):dim*i, 1+(dim*(j-1)):dim*j) = auxB;

        end
    end

    controlPoints = B'*W*B\B'*W*points;
    controlPoints = reshape(controlPoints, dim, []);
end

function b = binomial(x, y)
    b = factorial(x)/(factorial(y)*factorial(x-y));
end

function b = bernsteinPol(n, k, t)
    b = binomial(n, k)*t^k*(1-t)^(n-k);
end

function w = weights(k, L, times)
    w = zeros(L, 1);
    tL = times(L);
    for j = 1:L-1
        w(j) = tanh(k/(tL - times(j)));
    end
    w(end) = 1;
end

function [H, D] = buildHD(N, r, dim)
    H = zeros(N-r+1, N-r+1);
    
    for i=0:N-r
        for j = 0:N-r
            H(i+1, j+1) = (binomial(N-r, i)*binomial(N-r, j)/binomial(2*N-r, i+j)) / (2*N-r+1);
        end
    end
    H = kron(H, eye(dim));
    
    D = zeros(N-r+1, N+1);
    for i = 1:N-r+1
        D(i,i) = 1;
        D(i,i+1) = -2;
        D(i, (i+r)) = 1;
    end
    D = kron(D, eye(dim));
end

function B = buildB(n, tau, dim, L, w)
    B = 0;
    for i = 1:L
        b = zeros(1, n+1);
        for k = 0:n
            b(k+1) = bernsteinPol(n, k, tau(i));
        end
        b = kron(b, eye(dim));
        B = B + w(i)*b'*b;
    end
end

function Bp = buildBp(n, points, tau, dim, L, w)
    Bp = 0;
    for i = 1:L
        b = zeros(1, n+1);
        for k = 0:n
            b(k+1) = bernsteinPol(n, k, tau(i));
        end
        b = kron(b, eye(dim));
        Bp = Bp - w(i)*points(:, i)'*b;
    end
end