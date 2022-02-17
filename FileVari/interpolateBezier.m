clear all
clc


%% Generic 5th order Polynomial in the plane
% points = [1 3 6 8 10 20; 2 10 10 -11 12 10; 0.0 0.2 0.4 0.5 0.6 0.7];
% 
% controlPoints = interpolateBezierWLS(points, 5, eye(12));
% 
% t = 0:0.01:1;
% bezierCurve = zeros(2, length(t));
% for tt = 1:length(t)
%     bezierCurve(:, tt) = [0; 0];
%     for k = 1:6
%         bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(5, k-1, t(tt))*controlPoints(:, k);
%     end
% end
% 
% 
% figure(1)
% plot(bezierCurve(1, :), bezierCurve(2, :));
% hold on 
% plot(points(1, :), points(2, :), 'x');
% plot(controlPoints(1, :), controlPoints(2, :), 'o');
% hold off

%% Parabola in the Space
% tFinal = 3;
% t = linspace(0, tFinal, 100);
% 
% parabola = zeros(4, length(t));
% 
% parabola(1, :) = 3 + 0.5*t;
% parabola(2, :) = 3 + 0.5*t;
% parabola(3, :) = 2 + 15*t - 4.9*t.^2;
% parabola(4, :) = t/tFinal;
% 
% controlPoints = interpolateBezierLS(parabola(:, 1:5:30), 5);
% 
% t = 0:0.01:1;
% bezierCurve = zeros(3, length(t));
% for tt = 1:length(t)
%     bezierCurve(:, tt) = [0; 0; 0];
%     for k = 1:6
%         bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(5, k-1, t(tt))*controlPoints(:, k);
%     end
% end
% 
% figure(2)
% plot3(parabola(1, :), parabola(2, :), parabola(3, :));
% hold on
% plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))

%% Parabola with Noise
% rng(1)
% 
% dim = 3;
% tFinal = 3;
% t = linspace(0, tFinal, 8);
% 
% parabola = zeros(dim+1, length(t));
% 
% sigma = 0.05;
% sigmaMatrix = eye(dim*length(t)) * sigma^2;
% 
% parabola(1, :) = 3 + 0.5*t + sigma*randn(1, length(t));
% parabola(2, :) = 3 + 0.5*t + sigma*randn(1, length(t));
% parabola(3, :) = 2 + 15*t - 4.9*t.^2 + sigma*randn(1, length(t));
% parabola(4, :) = t/tFinal;
% 
% controlPoints1 = interpolateBezierLS(parabola, 5);
% controlPoints2 = interpolateBezierWLS(parabola, 5, sigmaMatrix);
% 
% t = 0:0.01:1;
% bezierCurve = zeros(3, length(t));
% for tt = 1:length(t)
%     bezierCurve(:, tt) = [0; 0; 0];
%     for k = 1:6
%         bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(5, k-1, t(tt))*controlPoints2(:, k);
%     end
% end
% 
% figure(3)
% plot3(parabola(1, :), parabola(2, :), parabola(3, :));
% hold on
% plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))

%% MONTECARLO PARABOLIC MOTION
m = 3;                  % dimension of the space
n = 5;                  % degree of Bezier Curve
tFinal = 10.0;          % time window length for motion estimation
tLast = 1.0;            % time window length for points collection
collectionFreq = 20;    % frequency of collection of the points
rng(1)                  % initialize the random number generator
sigma = 0.001;          % standard deviation of gaussian noise

% Initial Condition of the moving obstacle
p0 = [0; 0; 0];
v0 = [10; 10; 15];

% Collect Noisy Points
collectedPoints = zeros(m+1, cast(tLast*collectionFreq+1, 'uint8'));
collectedPoints(m+1, :) = 0:1/collectionFreq:tLast;

collectedPoints(1, :) = p0(1) + v0(1)*collectedPoints(4, :) + sigma*randn(1, length(collectedPoints(4, :))); 
collectedPoints(2, :) = p0(2) + v0(2)*collectedPoints(4, :) + sigma*randn(1, length(collectedPoints(4, :))); 
collectedPoints(3, :) = p0(3) + v0(3)*collectedPoints(4, :) - 4.9*collectedPoints(4, :).^2 + sigma*randn(1, length(collectedPoints(4, :)));

% Predict future points
obstacleDynamics = c2d(ss([zeros(m) eye(m); zeros(m) zeros(m)], [zeros(m+m-1, 1); 1], [eye(m,m), zeros(m,m)], []), 1/collectionFreq);

x = [collectedPoints(1:end-1, :); zeros(3,size(collectedPoints, 2))];
x(:, 1) = [collectedPoints(1:end-1, 1); zeros(3,1)];

for i=2:size(collectedPoints, 2)-1
    x(:, i+1) = obstacleDynamics.A*x(:, i) - [eye(m); eye(m)]*(obstacleDynamics.C*x(:, i) - collectedPoints(1:end-1, i+1)) - obstacleDynamics.B*9.8;
end

y = x(:, end);
for k = 1:30
    y(:, k+1) = obstacleDynamics.A*y(:, k) - obstacleDynamics.B*9.8;
end

y = y(1:3, 2:end);
y = [y; collectedPoints(end, end) + (1:1:size(y, 2))/collectionFreq];


collectedPoints(m+1, :) = collectedPoints(m+1, :)/tFinal;
y(end, :) = y(end, :) / tFinal;

% Fit Bezier Curve
controlPoints = interpolateBezierLS([collectedPoints, y], n);



t = 0:0.01:1;
bezierCurve = zeros(m, length(t));
for tt = 1:length(t)
    bezierCurve(:, tt) = zeros(m, 1);
    for k = 1:n+1
        bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(n, k-1, t(tt))*controlPoints(:, k);
    end
end

parabola = [ p0(1) + v0(1)*t*tFinal; ...
             p0(2) + v0(2)*t*tFinal; ...
             p0(3) + v0(3)*t*tFinal - 4.9*(t*tFinal).^2];

error = vecnorm(bezierCurve - parabola);        
figure(1)
plot3(parabola(1, :), parabola(2, :), parabola(3, :));
figure(2)
plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))
figure(3)
plot(error);

figure(4)
plot3(parabola(1, :), parabola(2, :), parabola(3, :));
hold on
plot3(y(1, :), y(2, :), y(3, :))

%plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))

legend('Parabola', 'Predicted Points', 'Bezier')

hold off
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
