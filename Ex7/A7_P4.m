clc
clear all
%% Assignment 7 Problem 4 %%

% Parameters
k1 = 1;
k2 = 1;
k3 = 1;
T = 0.1;
n = 2;
x0 = [5 1]';
xhat0 = [6 0]';
nx = 2;
nu = 1;
NN = 50;



% Continuous
Ac = [0 1;
    -k2, -k3];
Bc = [0 k3]';

% Discrete
Ad = eye(n) + T*Ac;
Bd = T*Bc;
C = [1 0];

% Optimization weights
Q = [4 0;0 4];
q = 4;
R = 1;

[K P e] = dlqr(Ad, Bd, Q/2, R/2, []);

%% MPC
% We want to use MPC control with time horizon of 10 steps-
% We need to generate equality constraints for 10 time steps

N = 10;

Q1 = q*ones(1,nx*(N-1));
Q1 = diag(Q1);
R1 = diag(R*ones(1,nu*N));
H = blkdiag(Q1,P,R1);


A1 = eye(nx*N);

for i = 1:nx:(N-1)*nx
    A1(i+2:i+3,i:i+1) = -Ad;
end

A2 =  kron(eye(N),-Bd);

Aeq = [A1 A2];

beq = zeros(N*nx,1);

xu = inf*ones(1,nx*N);
xl = -inf*ones(1,nx*N);

uu = 4*ones(1,nu*N);
ul = -4*ones(1,nu*N);

ub = [xu uu];
lb = [xl ul];

% Allocate space for solutions
x = zeros(nx,NN+1);
xhat = zeros(nx,NN+1);
u = zeros(nu,NN+1);

%Initial value

x(:,1) = x0;

for i = 1:NN
    % Update initial value
    beq(1:nx,1) = Ad*x(:,i);
    
    % Solve finite time horizon optimization
    [z,fval,exitflag,output,lambda] = quadprog(H,[],[],[],Aeq,beq,lb,ub);
    
    % Get input value
    u(i) = z(N*nx+1);
    
    %Simulate one additional step
    x(:,i+1) = Ad*x(:,i) + Bd*u(i);
    
end
    


% Plotting solution

t0 = 0;
t1 = NN*T;

t = t0:T:t1;

figure(2)

xlabel('t [s]');
subplot(211)
plot(t,x,'-');
hleg1 = legend('$x_1(t)$','$x_2(t)$');
set(hleg1, 'Interpreter', 'Latex');
grid on

subplot(212)
plot(t,u)
hleg2 = legend('$u(t)$');
set(hleg2, 'Interpreter', 'Latex')
grid on

