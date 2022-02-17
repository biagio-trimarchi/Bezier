% REPLANNING BEZIER
% Notes

%% Clear workspace
clear 
clc 
addpath('./Functions')

%% PARAMETERS
d = 2;
n = 5;
dT = 2;
N = 2;
Tf = N*dT;
t = 0:0.01:1;   % Time Vector
dSafe = 1;

%% COMPUTE Q
Q = computeQ(n);

%% DRONES CONTROL POINTS
P1 = [0,0,0,1,2,3; 
     0,0,0,3,1,4;
     3,4,5,  6,7,7;
     4,7,15, 7,4,4];

P2 = [3, 3, 3, 2, 1, 3.5;  
      0, 0, 0, 1, 2, 4.5;
      3.5, 5.5, 11, 3, 2, 4;
      4.5, 6.5, 7, 2, 1, 0] ;

  
%% OPTIMIZATION PROBLEM
% Equality Constraints
Aaux1 = [1 0 0 0 0 0;
        -1 1 0 0 0 0;
        1 -2 1 0 0 0];
    
Aaux2 = [0 0 0 0 0 1;
        0 0 0 0 -1 1;
        0 0 0 1 -2 1];

Aaux3 = [1 0 0 0 0 0;
        -1 1 0 0 0 0;
         1 -2 1 0 0 0;
         0 0 0 0 0 1;
         0 0 0 0 -1 1;
         0 0 0 1 -2 1];
     
Aeq = [blkdiag(Aaux1, Aaux2, Aaux1, Aaux2); blkdiag([Aaux2, -Aaux1], [Aaux2, -Aaux1])];
%Aeq = [Aeq; -Aaux2, Aaux1; -Aaux2, Aaux1];

beq = [P1(1, 1);
       P1(1, 2) - P1(1, 1);
       P1(1, 3) - 2*P1(1, 2) + P1(1,1);
       P1(3, end);
       P1(3, end) - P1(3, end-1);
       P1(3, end-2) - 2*P1(3, end-1) + P1(3,end);
       P1(2, 1);
       P1(2, 2) - P1(2, 1);
       P1(2, 3) - 2*P1(2, 2) + P1(2,1);
       P1(4, end);
       P1(4, end) - P1(4, end-1);
       P1(4, end-2) - 2*P1(4, end-1) + P1(4,end);
       0;
       0;
       0;
       0;
       0;
       0;];

% Dynamic Constraints ...
for i =1:n
    [zeros(1, i-1), -1, 1, zeros(1, n+1-2-(i-1))];
    AdynV(i, :) = [zeros(1, i-1), -1, 1, zeros(1, n+1-2-(i-1))];
end

for i =1:n-1
    AdynA(i, :) = [zeros(1, i-1), 1, -2, 1, zeros(1, n+1-3-(i-1))];
end

Adyn = [blkdiag(AdynV, AdynV, AdynV, AdynV); blkdiag(AdynA, AdynA, AdynA, AdynA)];
Adyn = [Adyn; -Adyn];

bDynV = 20*ones(n, 1);
bDynA = 20*ones(n-1, 1);
bDyn = [bDynV; bDynV; bDynV; bDynV; bDynA; bDynA; bDynA; bDynA];
bDyn = [bDyn; bDyn];

% Functions
J = @(x) costP(x, n, d, P1(1:2, n+1));
constr = @(x) nlCnstP(x, P2(1:2, :), n, d, dSafe, Q);

% Initial guess
x0 = [P1(1,:), P1(3, :), P1(2, :), P1(4, :)];

opts = optimoptions('fmincon','Algorithm','sqp');
tic
[x, fVal] = fmincon(J, x0, Adyn, bDyn, Aeq, beq, [], [], constr, opts);
toc

P3 = [x(1:n+1); x((n+1)*d+1:(n+1)*(d+1)); x(n+2:d*(n+1)); x((d+1)*(n+1)+1:end)];

%% PLOT CURVES

bezierCurve11 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve11(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve11(:, tt) = bezierCurve11(:, tt) + bernsteinPol(n, k-1, t(tt))*P1(1:2, k);
    end
end

bezierCurve12 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve12(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve12(:, tt) = bezierCurve12(:, tt) + bernsteinPol(n, k-1, t(tt))*P1(3:4, k);
    end
end

bezierCurve21 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve21(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve21(:, tt) = bezierCurve21(:, tt) + bernsteinPol(n, k-1, t(tt))*P2(1:2, k);
    end
end

bezierCurve22 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve22(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve22(:, tt) = bezierCurve22(:, tt) + bernsteinPol(n, k-1, t(tt))*P2(3:4, k);
    end
end

bezierCurve31 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve31(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve31(:, tt) = bezierCurve31(:, tt) + bernsteinPol(n, k-1, t(tt))*P3(1:2, k);
    end
end

bezierCurve32 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve32(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve32(:, tt) = bezierCurve32(:, tt) + bernsteinPol(n, k-1, t(tt))*P3(3:4, k);
    end
end

bezierCurve12 = zeros(d, length(t));
for tt = 1:length(t)
    bezierCurve12(:, tt) = zeros(d, 1);
    for k = 1:n+1
        bezierCurve12(:, tt) = bezierCurve12(:, tt) + bernsteinPol(n, k-1, t(tt))*P1(3:4, k);
    end
end


figure(1)
plot(bezierCurve11(1, :), bezierCurve11(2, :));
hold on
plot(bezierCurve12(1, :), bezierCurve12(2, :));
plot(bezierCurve21(1, :), bezierCurve21(2, :));
plot(bezierCurve22(1, :), bezierCurve22(2, :));
plot(P2(1, :), P2(2, :), 'o');
hold off

figure(2)
plot(bezierCurve31(1, :), bezierCurve31(2, :));
hold on
plot(bezierCurve32(1, :), bezierCurve32(2, :));
plot(bezierCurve21(1, :), bezierCurve21(2, :));
plot(bezierCurve22(1, :), bezierCurve22(2, :));
plot(P2(1, :), P2(2, :), 'o');
hold off


for i = 1:2
    figure(11)
    if i == 1
        for tt =1:length(t)
            
            plot(bezierCurve11(1, :), bezierCurve11(2, :));
            hold on
            plot(bezierCurve21(1, :), bezierCurve21(2, :));
            plot(bezierCurve12(1, :), bezierCurve12(2, :));
            plot(bezierCurve22(1, :), bezierCurve22(2, :));
            
            %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
            plotCircle(bezierCurve11(1, tt), bezierCurve11(2, tt), dSafe)
            plot(bezierCurve11(1, tt), bezierCurve11(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
            plot(bezierCurve21(1, tt), bezierCurve21(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
            hold off
            
            xlim([-1, 11])
            ylim([-1, 11])
            title("Simulation New Trajectory")
            drawnow
        end
    else
        for tt =1:length(t)
            plot(bezierCurve11(1, :), bezierCurve11(2, :));
            hold on
            plot(bezierCurve21(1, :), bezierCurve21(2, :));
            plot(bezierCurve12(1, :), bezierCurve12(2, :));
            plot(bezierCurve22(1, :), bezierCurve22(2, :));
            
            
            %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
            plotCircle(bezierCurve12(1, tt), bezierCurve12(2, tt), dSafe)
            plot(bezierCurve12(1, tt), bezierCurve12(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
            plot(bezierCurve22(1, tt), bezierCurve22(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
            hold off
            
            xlim([-1, 11])
            ylim([-1, 11])
            title("Simulation New Trajectory")
            drawnow
        end
    end
end

for i=1:2
    figure(12)
    if i == 1
        for tt =1:length(t)
            
            plot(bezierCurve31(1, :), bezierCurve31(2, :));
            hold on
            plot(bezierCurve21(1, :), bezierCurve21(2, :));
            plot(bezierCurve32(1, :), bezierCurve32(2, :));
            plot(bezierCurve22(1, :), bezierCurve22(2, :));
            
            %plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
            %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
            %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
            plotCircle(bezierCurve31(1, tt), bezierCurve31(2, tt), dSafe)
            plot(bezierCurve31(1, tt), bezierCurve31(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
            plot(bezierCurve21(1, tt), bezierCurve21(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
            hold off
            
            xlim([-1, 11])
            ylim([-1, 11])
            title("Simulation New Trajectory")
            drawnow
        end
    else
        for tt =1:length(t)
            plot(bezierCurve31(1, :), bezierCurve31(2, :));
            hold on
            plot(bezierCurve21(1, :), bezierCurve21(2, :));
            plot(bezierCurve32(1, :), bezierCurve32(2, :));
            plot(bezierCurve22(1, :), bezierCurve22(2, :));
            
            %plot(bezierCurved(1, tt), bezierCurved(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
            %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), safeDist)
            %plotCircle(bezierCurved(1, tt), bezierCurved(2, tt), rDetect)
            plotCircle(bezierCurve32(1, tt), bezierCurve32(2, tt), dSafe)
            plot(bezierCurve32(1, tt), bezierCurve32(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0]);
            plot(bezierCurve22(1, tt), bezierCurve22(2, tt), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1]);
            hold off
            
            xlim([-1, 11])
            ylim([-1, 11])
            title("Simulation New Trajectory")
            drawnow
        end
    end
end