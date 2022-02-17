clear 
clc
addpath('./Functions')

%% Parameters
n = 5;          % Spatial Curve Order
m = 3;          % Temporal Parametrization Order
d = 2;          % Space Order
t = 0:0.01:1;   % Time Vector

safeDist = 0.25;    % Safe Distance 
vMax = 10;          % Maximum Velocity
aMax = 5;          % Maximum Acceleration
tF = 50;            % Total Time

%% Divide Curve
Pd(:, 1) = [0.0; 0.0];
Pd(:, 2) = [0.2; 1.0];
Pd(:, 3) = [0.5; 2.4];
Pd(:, 4) = [1.2; 1.1];
Pd(:, 5) = [1.5; 2.0];
Pd(:, 6) = [2.0; 3.0];

bezierCurve = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(n, k-1, t(tt))*Pd(:, k);
    end
end

[P1x, P2x, Pdivx] = bezierDivision(Pd(1, :), 0.8, n);
[P1y, P2y, Pdivy] = bezierDivision(Pd(2, :), 0.8, n);

P1 = [P1x'; P1y'];
P2 = [P2x'; P2y'];

%% Plots

plot(Pd(1, :), Pd(2,:), 'o')
hold on 
plot(P1(1, :), P1(2,:), 'x')
plot(P2(1, :), P2(2,:), 's')
plot(bezierCurve(1, :), bezierCurve(2, :))
hold off

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

% Divide Curve