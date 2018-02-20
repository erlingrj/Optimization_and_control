clc
clear all

%Problem 2b Model Predictive Control. From Excercise 5

r = 1;
N= 30;
x0 = [0 0 1]';
nx = 3;
nu = 1;

b_length = ([1 1 2 4 8 14]);
n = length(b_length);

% System matrices
Ad = [0 0 0;0 0 1;0.1 -0.79 1.78];
Bd = [1 0 0.1]';
Cd = [0 0 1];

%Cost function weights
Q = [0 0 0;0 0 0;0 0 1];
R = r;

%Objective function. 0.5z'Gz
G1 = kron(eye(N),2*Q);
G2 = kron(eye(n),2*R);
G = blkdiag(G1,G2);

% Equality constraint
Aeq_c1 = eye(N*nx);                             % Component 1 of A_eq
Aeq_c2 = kron(diag(ones(N-1,1),-1), -Ad);        % Component 2 of A_eq
ones_block = blkdiag(ones(b_length(1),1), ...
                     ones(b_length(2),1), ...
                     ones(b_length(3),1), ...
                     ones(b_length(4),1), ...
                     ones(b_length(5),1), ...
                     ones(b_length(6),1));      % Block-diagonal matrix of 1-vectors
Aeq_c3 = kron(ones_block, -Bd);                  % Component 3 of A_eq
Aeq = [Aeq_c1 + Aeq_c2, Aeq_c3];

%Upper and lower bounds
x_lb = -Inf(3*N,1);
x_ub = Inf(3*N,1);

u_lb = -1*ones(n,1);
u_ub = ones(n,1);

lb = [x_lb; u_lb];
ub = [x_ub; u_ub];

x = zeros(nx,N+1);
u = zeros(nu,N);
x(:,1) = x0;


beq = zeros(nx*N,1);

for i = 1:N
    beq(1:3,1) = Ad*x(:,i);
    
    [sol,fval,exitflag,output] = quadprog(G,[],[],[],Aeq,beq,lb,ub);
    
    x(:,i) = sol(1:3);
    u_blocks = sol(nx*N+1:nx*N+n);

    u(i) = u_blocks(1);
    
    x(:,i+1) = Ad*x(:,i) + Bd*u(i);

end