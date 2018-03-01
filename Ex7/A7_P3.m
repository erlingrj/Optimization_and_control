clc
clear all
% Exercise 7 Problem 3 %

%% MPC control %%
%Parameters
k1 = 1;
k2 = 1;
k3 = 1;
T = 0.1;
n = 2;
x0 = [5 1]';
xhat0 = [6 0]';
nx = 2;
nu = 1;



%Contineous
Ac = [0 1;
    -k2, -k3];
Bc = [0 k3]';

% Discrete
Ad = eye(n) + T*Ac;
Bd = T*Bc;
C = [1 0];

% Optimization weights
q = 4;
R = 1;


% We want to use MPC control with time horizon of 10 steps-
% We need to generate equality constraints for 10 time steps

N = 10;

Q1 = q*ones(1,nx*N);
Q1 = diag(Q1);
R1 = diag(R*ones(1,nu*N));
H = blkdiag(Q1,R1);


A1 = eye(nx*N);

for i = 1:nx:(N-1)*nx
    A1(i+2:i+3,i:i+1) = -Ad;
end

A2 =  kron(eye(N),-Bd);

Aeq = [A1 A2];

beq = zeros(N*nx,1);
beq(1:nx,1) = Ad*x0;

xu = inf*ones(1,nx*N);
xl = -inf*ones(1,nx*N);

uu = 4*ones(1,nu*N);
ul = -4*ones(1,nu*N);

ub = [xu uu];
lb = [xl ul];

[z,fval,exitflag,output,lambda] = quadprog(H,[],[],[],Aeq,beq,lb,ub);

