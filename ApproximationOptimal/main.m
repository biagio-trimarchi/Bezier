% OPTIMAL CONTROL BEZIER APPROXIMATION EXAMPLE
% notes

%% Clear Workspace
clear 
clc 
addpath('./Functions')

%% Parameters
d = 2;          % Space Dimension
N = 5;          % Bezier Order Approximation
E = 50;         % Safe Distance
Cp = 5*d/4;
dP = Cp/N;      % Relaxation
w = 1/(N+1);    % Parameter
t = 0:0.01:1;   % Time Vector
%% Obstacles
pO = [   0,  450,  850; ...
      -800, -750, -730];

%% Input Limits
uMin = 15;
uMax = 32;

%% Initial And Final Conditions
x0 = [-500; -900];
xF = [1500; -600];

%% Optimization Problem
nlc = @(xut) nlCons(xut, N, d, uMax, uMin, E, pO, 3, 0);

J = @(x) x(end);

Aeq(1:d, :) = [eye(d), zeros(d, 2*N+2*(N+1)+1)];
Aeq(d+1:d+d, :) = [zeros(d, 2*N), eye(d), zeros(d, 2*(N+1)+1)];

beq(1:d) = x0;
beq(d+1:d+d) = xF;

lb = [-inf*ones(1, 4*(N+1)), 0];

opts = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e4);
xut = fmincon(J, zeros(4*(N+1)+1,1), [], [], Aeq, beq, lb, [], nlc, opts);

%% Plots

x = reshape(xut(1:(N+1)*d), d, []);
u = reshape(xut((N+1)*d+1:end-1), d, []);
tF = xut(end)*w;

bezX = zeros(d, length(t));

for tt=1:length(t)
    bezX(:, tt) = bezier(N, t(tt), x, d);
end


bezU = zeros(d, length(t));

for tt=1:length(t)
    bezU(:, tt) = bezier(N, t(tt), u, d);
end

figure(1)
plot(bezX(1, :), bezX(2,:))
hold on
plot(x(1, :), x(2,:), 'x')
plot(pO(1, :), pO(2, :), 'o')
plotCircle(pO(1, 3), pO(2, 3), E);
hold off

