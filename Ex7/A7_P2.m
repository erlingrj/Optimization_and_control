% Excercise 7
%% Problem 2
clc
clear all


%Parameters
k1 = 1;
k2 = 1;
k3 = 1;
T = 0.1;
n = 2;
x0 = [5 1]';
xhat0 = [6 0]';


%Contineous
Ac = [0 1;
    -k2, -k3];
Bc = [0 k3]';

% Discrete
Ad = eye(n) + T*Ac;
Bd = T*Bc;
C = [1 0];

%% Optimization weights
Q = [4 0;0 4];
R = 1;
N = [];

% Linear-quadratic state-feedback solver

[K S e] = dlqr(Ad,Bd,Q,R,N);

% Eigenvalues and eigenvectors of the closed loop system.
[V1 D1] = eig(Ad-Bd*K);

% Placing poles for state observer.
p = [0.5+0.03i; 0.5-0.03i];
L = place(Ad',C',p)';

% Simulate for N = 50
N = 50;
x = zeros(2,N+1);
xhat = zeros(2,N+1);
u = zeros(1,N+1);
y = zeros(1,N+1);
yhat = zeros(1,N+1);

% Initial values
x(:,1) = x0;
xhat(:,1) = xhat0;
y(1) = C*x(:,1);
yhat(1) = C*xhat(:,1);

for i = 1:N-1
    % Calculate optimal gain
    u(i) = -K*xhat(:,i);
    
    % Simulate one step ahead
    x(:,i+1) = Ad*x(:,i) + Bd*u(i);
    y(i+1) = C*x(:,i+1);
    
    
    %Update estimate
    xhat(:,i+1) = Ad*xhat(:,i) + L*(y(i) - yhat(i)) + Bd*u(i);
    yhat(i+1) = C*xhat(:,i+1);
end

t0 = 0;
t1 = N*T;

t = 0:T:t1;

figure(1)
hold on
xlabel('t [s]')
plot(t,x, '-')
plot(t,xhat,'--')
hleg = legend('$x_1(t)$','$x_2(t)$','$\hat{x}_1(t)$','$\hat{x}_2(t)$');
set(hleg, 'Interpreter', 'Latex')


    



