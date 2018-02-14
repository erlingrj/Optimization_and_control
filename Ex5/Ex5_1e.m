clc
clear all


r = 1;
N= 30;
x0 = [0 0 1]';


% System matrices
Ad = [0 0 0;0 0 1;0.1 -0.79 1.78];
Bd = [1 0 0.1]';
Cd = [0 0 1];

%Cost function weights
Q = [0 0 0;0 0 0;0 0 2];
R = r;

%Objective function. 0.5z'Gz
G1 = kron(eye(N),Q);
G2 = kron(eye(N),R);
G = blkdiag(G1,G2);

f = [];

% Equality constraints Aeq = beq
A1 = eye(3*N);

for i = 3:3:(3*(N-1));
    A1(i+1:i+3,i-2:i) = -Ad;
end
     
A2 =  kron(eye(N),-Bd);
Aeq = [A1 A2];

beq = zeros(3*N,1);
beq(1:3,1) = Ad*x0;

% Solve quadratic optimization problem
sol = quadprog(G,f,[],[],Aeq,beq);

z = sol(1:(4*N));
u = z(3*N + 1:4*N);
y = [x0(3); z(3:3:3*N)];

t = 1:N;

figure
subplot(211)
plot([0 t],y)
legend y(t)
xlabel('t[s]')
ylabel('y(t)')

subplot(212)
plot(t,u)
legend u(t)
xlabel('t[s]')
ylabel('u(t)')


