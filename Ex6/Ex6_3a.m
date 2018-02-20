clc
clear all

% Plant model
A = [0 0 0;
    0 0 1;
    0.1 -0.79 1.78];
B = [1 0 0.1]';
C = [0 0 1];
x0 = [0 0 1]';

mx = 3;
mu = 1;
N = 30;

% Cost/objective function
Q = diag([0 0 1]);
R = 1;

H = blkdiag(kron(eye(N),Q),kron(eye(N),R));

% Equality constraints
Aeq1 = eye(mx*N);

for i = 3:3:((N-1)*mx)
    Aeq1(i+1:i+3,i-2:i) = -A;
end

Aeq2 = kron(eye(N*mu),-B);
Aeq = [Aeq1 Aeq2];

beq = zeros(1,mx*N);
beq(1,1:mx) = A*x0;

% Upper and lower bounds
xu = inf*ones(N*mx, 1);
xl = -inf*ones(N*mx, 1);

uu = ones(N*mu, 1);
ul = -1*ones(N*mu, 1);

vu = [xu ; uu];
vl = [xl ; ul];

% Solve LQ problem
[X,FVAL,EXITFLAG] = quadprog(H,[],[],[], Aeq, beq, vl, vu);

%Plot solution
z = X(1:(4*N));
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
