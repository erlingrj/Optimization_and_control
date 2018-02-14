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

% Equality constraints Aeq = beq
A1 = eye(3*N);

for i = 3:3:(3*(N-1));
    A1(i+1:i+3,i-2:i) = -Ad;
end
     
A2 =  kron(eye(N),-Bd);
Aeq = [A1 A2];

beq = zeros(3*N,1);
beq(1:3,1) = Ad*x0;

% KKT System.

zero_matrix = zeros(3*N);
zero_vector = zeros(4*N,1);

KKT_matrix = [G -Aeq';Aeq zero_matrix];
KKT_vector = [zero_vector;beq];

KKT_sol = KKT_matrix\KKT_vector;

z = KKT_sol(1:(4*N));
u = z(3*N + 1:4*N);
y = [x0(3); z(3:3:3*N)];

t = 0:N;

figure
plot(t,y)
legend y(t)
xlabel('t[s]')
ylabel('y(t)')