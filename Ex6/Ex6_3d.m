clc
clear all

%Problem 2b Model Predictive Control. From Excercise 5

r = 1;
N= 30;
x0 = [0 0 1]';


% System matrices
Ad = [0 0 0;0 0 1;0.1 -0.79 1.78];
Bd = [1 0 0.1]';
Cd = [0 0 1];

%Cost function weights
Q = [0 0 0;0 0 0;0 0 1];
R = r;

%Objective function. 0.5z'Gz
G1 = kron(eye(N),2*Q);
G2 = kron(eye(N),2*R);
G = blkdiag(G1,G2);

% Equality constraints Aeq = beq
A1 = eye(3*N);

%Upper and lower bounds
x_lb = -Inf(3*N,1);
x_ub = Inf(3*N,1);

u_lb = -1*ones(N,1);
u_ub = ones(N,1);

lb = [x_lb; u_lb];
ub = [x_ub; u_ub];

for i = 3:3:(3*(N-1))
    A1(i+1:i+3,i-2:i) = -Ad;
end
     
A2 =  kron(eye(N),-Bd);
Aeq = [A1 A2];
beq = zeros(3*N,1);

x = zeros(3,N+1);
u = zeros(1,N);
x(:,1) = x0;

for i = 1:N
    beq(1:3,1) = Ad*x(:,i);
    
    [sol,fval,exitflag,output] = quadprog(G,[],[],[],Aeq,beq,lb,ub);
    
    x(:,i) = sol(1:3);
    u(i) = sol(3*N + 1);
    
    x(:,i+1) = Ad*x(:,i) + Bd*u(i);

end