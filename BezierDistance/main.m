% BEZIER DISTANCE 

%% Clear workspace
clear 
clc 
addpath('./Functions')

P1 = [0 1 3 
      0 1 2 
      0 0 0];
 
Q1 = [2 3 1
      3 3 1
      0 0 0];
[x, P2, P3] = castel(P1, 1/2);
 
tic
computeDistance(P1, Q1, 10, 0.5)
toc
