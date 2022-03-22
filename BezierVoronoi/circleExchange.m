%% %% %% BEZIER AND VORONOI
% AAAA

%% SETUP
clear
clc
addpath('./Functions')

%% PARAMETER
d = 2;  % Space Dimension
r = 5;  % Bezier Degree
N = 20;  % Number of Agents

%% AGENT INITIALIZATION

R = 20; % Circle Radius
for i = 1:N
    agent(i).x0 = [R*cos(i/N * 2*pi); R*sin(i/N * 2*pi)];
    agent(i).xF = -agent(i).x0;
end

figure(1)
plot(0, 0)
hold on
for i = 1:N
    plot(agent(i).x0(1), agent(i).x0(2), 'x', 'MarkerSize', 10)
end

%% RECEDING HORIZION PARAMETERS

nW = 5;     
T = 1;
