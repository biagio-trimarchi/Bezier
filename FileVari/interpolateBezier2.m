clear all
clc

%% TEST
L = 15;
degree = 5;
n = 3;

tFinal = 10;
t = 0:0.5:tFinal;

parabola = zeros(n+1, length(t));

parabola(1, :) = 3 + 0.5*t;
parabola(2, :) = 3 + 0.5*t;
parabola(3, :) = 2 + 15*t - 4.9*t.^2;
parabola(4, :) = t;

points = parabola(:, 1:L);

tL = points(n+1, end);

w = weights(1, L, points(4,:));
w = ones(L,1);
[H, D] = buildHD(degree, 2, n);

J2 = D'*H*D;

B = buildB(degree, points(4, :)/tFinal, n, L, w);
Bp = buildBp(degree, points(1:3, :), points(4, :)/tFinal, n, L, w);

% Interpolate 1
points(4, :) = points(4, :)/tFinal;
controlPoints1 = interpolateBezierLS(points, degree);
t = 0:0.01:1;
bezierCurve = zeros(3, length(t));
for tt = 1:length(t)
    bezierCurve(:, tt) = [0; 0; 0];
    for k = 1:6
        bezierCurve(:, tt) = bezierCurve(:, tt) + bernsteinPol(5, k-1, t(tt))*controlPoints1(:, k);
    end
end

% Interpolate 2
a = [1 -2 1; -1 2 -1];
A = kron(a, eye(degree+1));
b = [10/(degree*(degree-1));10/(degree*(degree-1))];
b = kron(b, ones(degree+1, 1));

[controlPoints2, fval] = quadprog(B + 0.5*J2, Bp, [], []);
controlPoints2 = reshape(controlPoints2, 3, []);


bezierCurve2 = zeros(3, length(t));
for tt = 1:length(t)
    bezierCurve2(:, tt) = [0; 0; 0];
    for k = 1:6
        bezierCurve2(:, tt) = bezierCurve2(:, tt) + bernsteinPol(5, k-1, t(tt))*controlPoints2(:, k);
    end
end

figure(1)
%plot3(parabola(1, :), parabola(2, :), parabola(3, :));
%hold on
plot3(bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :))
hold on
plot3(bezierCurve2(1, :), bezierCurve2(2, :), bezierCurve2(3, :))
%% QUADRATIC PROGRAM
% p = reshape(points(1:end-1, :), [], 1);
% b = zeros(n*(degree+1), 1);
% B = kron(eye(L), b*b');
% 
% for i = 1:L
%     for j = 0:degree
%         b(j*n+1:j*n+n) = bernsteinPol(degree, j, points(n+1, i)/tFinal)*ones(n, 1);
%     end
%     B(i:i+11, i:i+11) = b*b';
% end




%% EXTERNAL FUNCTIONS
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

function b = binomial(x, y)
    % Compute recursively the binomial (x y)
    
    b = factorial(x)/(factorial(y)*factorial(x-y));
end

function b = bernsteinPol(n, k, t)
    % Compute the value of the k-th Bernstein Polynomial of
    % of order n at time t
    
    b = binomial(n, k)*t^k*(1-t)^(n-k);
end

function b = bezier(n, t, controlPoints, dim)
    % Compute the coordinates of the Bezier curve of order n
    % at time t given the control points and the dimension of the space
    
    b = zeros(dim, 1);
    for k=0:n
        b = b + bernsteinPol(n, k, t)*controlPoints(:, k);
    end
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
        for j = 0:r-1
            D(i,i+j) = (-1)^(r+j)*binomial(r, j);
        end
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