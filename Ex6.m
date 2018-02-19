clc
clear all
% Task 1d Solve the stationary Riccati equation.

%System
Ad = [1 0.5; 0 1];
Bd = [0.125 0.5]';

% Cost function
Q = [2 0;0 2];
R = 2;

%Solve the Riccati function and fin optimal gain.
[K,S,e] = dlqr(Ad,Bd,Q/2,R/2); 