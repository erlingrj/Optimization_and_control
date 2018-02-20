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
n=6; %Number of "input-blocks"
input_length = [1 1 2 4 8 14];


% Cost/objective function
Qt = 2*diag([0 0 1]);
Rt = 2*1;

Q = kron(eye(N), Qt);
R = kron(N/n * eye(n), Rt);


H = blkdiag(Q,R);

% Equality constraints
Aeq1 = eye(mx*N);

for i = 3:3:((N-1)*mx)
    Aeq1(i+1:i+3,i-2:i) = -A;
end

% Change the B-part of Aeq to only allow a constant U over five steps

Aeq2 = zeros(6,N*mx);

offset = 1;
for i = 1:n
   for j = 1:input_length(i)
       Aeq2(i,offset:offset+mx-1) = -B;
       offset = offset + mx;
    end
end

Aeq = [Aeq1 Aeq2'];

beq = zeros(mx*N,1);
beq(1:mx,1) = A*x0;

% Upper and lower bounds
xu = inf*ones(N*mx, 1);
xl = -inf*ones(N*mx, 1);

uu = ones(n*mu, 1);
ul = -1*ones(n*mu, 1);

vu = [xu ; uu];
vl = [xl ; ul];

% Solve LQ problem
[y,fval,exitflag,output] = quadprog(H,[],[],[], Aeq, beq, vl, vu);

%Plot solution
z = y(1:length(y));
u = zeros(N,1);
offset = 1;
for i = 1:n
    u(offset:offset+input_length(i)-1) = z(mx*N+i)*ones(input_length(i),1);
    offset = offset + input_length(i);
end

y = [x0(3); z(3:3:3*N)];

t = 1:N;

figure
subplot(211)
plot([0 t],y,'-o')
legend y(t)
xlabel('t[s]')
ylabel('y(t)')

subplot(212)
plot(t,u,'-o')
legend u(t)
xlabel('t[s]')
ylabel('u(t)')
