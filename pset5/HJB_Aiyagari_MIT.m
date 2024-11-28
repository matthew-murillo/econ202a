%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: Aiyagari MIT Shock
% 
% Author: Matt Murillo
%
%
% Reference: Andreas Schaab, Kiyea Jin, aiyagari_poisson_MITshock.m by Benjamin Moll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();
dr = .0001;

% Storing excess demand for Newton's Method;
Z_r = [];
r_r = [];

% Productivity evolution
T = 200;
N = 400;
dt = T/N;
time = (0:N-1)*dt;

nu = .2;
A_t = zeros(N,1);
A_t(1)=.97*p.A;
for n=1:N-1
A_t(n+1) = dt*nu*([p.A-A_t(n)]) + A_t(n);
end

%% 2. INITIALIZE GRID POINTS

a = linspace(p.amin, p.amax, p.I)';
da = (p.amax-p.amin)/(p.I-1);

aa = [a, a]; % I*2 matrix

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward differential operator 
% Df such that Df*V=dVf and Db such that Db*V=dVb

Df = zeros(p.I, p.I);
for i = 1:p.I-1
    Df(i,i) = -1/da; Df(i,i+1) = 1/da;
end
Df = sparse(Df);

Db = zeros(p.I, p.I);
for i = 2:p.I
    Db(i,i-1) = -1/da; Db(i,i) = 1/da;
end
Db = sparse(Db);

% 3-2. Construct A_switch matrix

A_switch = [speye(p.I).*(-p.lambda(1)), speye(p.I).*p.lambda(1);
speye(p.I).*p.lambda(2), speye(p.I).*(-p.lambda(2))];

% 3-3. Setting rent, wage, labor, and capital levels
r = 0.043; % Initial guess
r_r = [r_r;r];
for i = 1:p.maxit
% for i = 1:10
p.r = r;
L = (p.z_u*p.lambda_e + p.z_e*p.lambda_u)/(sum(p.lambda));
K = (p.alpha*p.A/(p.r + p.d))^(1/(1-p.alpha))*L; % maybe problems here
w = (1-p.alpha)*p.A*K.^p.alpha*L^(-p.alpha);

% 3-4. Guess an initial value of the value function
zz = ones(p.I, 1).*p.zz; % I*2 matrix

% The value function of "staying put" 
v0 = p.u(w*zz + max(r,0.01).*aa)./p.rho; 
V = v0;


%% 4. VALUE FUNCTION ITERATION
    
    
for n=1:p.maxit
    
    % 4-1. Compute the derivative of the value function 
    dVf = Df*V;
    dVb = Db*V;
    
    % 4-2. Boundary conditions
    dVb(1,:) = p.mu(w*zz(1,:) + p.r.*aa(1,:)); % a>=a_min is enforced (borrowing constraint)
    dVf(end,:) = p.mu(w*zz(end,:) + p.r.*aa(end,:)); % a<=a_max is enforced which helps stability of the algorithm
    
    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)
    
    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf);
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = w*zz + p.r.*aa - cf;
    sb = w*zz + p.r.*aa - cb;
    
    % 4-5. Upwind scheme
    If = sf>0;
    Ib = sb<0;
    I0 = 1-If-Ib;
    dV0 = p.mu(w*zz + p.r.*aa); % If sf<=0<=sb, set s=0
    
    dV_upwind = If.*dVf + Ib.*dVb + I0.*dV0;
    
    c = p.inv_mu(dV_upwind);
    
    % 4-6. Update value function: 
    % Vj^(n+1) = [(rho + 1/Delta)*I - (Sj^n*Dj^n+A_switch)]^(-1)*[u(cj^n) + 1/Delta*Vj^n]
    
    V_stacked = V(:); % 2I*1 matrix
    c_stacked = c(:); % 2I*1 matrix
    
    % A = SD
    SD_u = spdiags(If(:,1).*sf(:,1), 0, p.I, p.I)*Df + spdiags(Ib(:,1).*sb(:,1), 0, p.I, p.I)*Db; % I*I matrix
    SD_e = spdiags(If(:,2).*sf(:,2), 0, p.I, p.I)*Df + spdiags(Ib(:,2).*sb(:,2), 0, p.I, p.I)*Db; % I*I matrix
    SD = [SD_u, sparse(p.I, p.I);
    sparse(p.I, p.I), SD_e]; % 2I*2I matrix
    
    % P = A + A_switch
    P = SD + A_switch;
    
    % B = [(rho + 1/Delta)*I - P]
    B = (p.rho + 1/p.Delta)*speye(2*p.I) - P; 
    
    % b = u(c) + 1/Delta*V
    b = p.u(c_stacked) + (1/p.Delta)*V_stacked;
    
    % V = B\b;
    V_update = B\b; % 2I*1 matrix
    V_change = V_update - V_stacked;
    V = reshape(V_update, p.I, 2); % I*2 matrix
    % V = .9*V + .1*reshape(V_update, p.I, 2); % I*2 matrix
    
    
    % 3-6. Convergence criterion
    dist(n) = max(abs(V_change));
    if dist(n)<p.tol
        disp('Value function converged. Iteration = ')
        disp(n)
        break
    end
end

toc;

%% 5. KF EQUATION



PT = P';
b = zeros(2*p.I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*p.I-i_fix)];
PT(i_fix,:) = row;

%Solve linear system
gg = PT\b;
g_sum = gg'*ones(2*p.I,1)*da;
gg = gg./g_sum;

g = [gg(1:p.I),gg(p.I+1:2*p.I)];

check1 = g(:,1)'*ones(p.I,1)*da;
check2 = g(:,2)'*ones(p.I,1)*da;

sum(check1 + check2)

gg = g;

%% Steady state employment and unemployment

T = p.lambda;
T(2) = -T(2);
T = [1,0; T];
b = [0.1;0];
s_hat = T\b;
s = s_hat./(sum(s_hat));

%% Asset supply
S = gg(:,1)'*a*da + gg(:,2)'*a*da;

%% Excess demand. Newton's Method
Z = K - S;
Z_r = [Z_r; Z];


if abs(Z) < p.tol
    break
end

if i == 1
    r = r + dr;
    r_r = [r_r;r]
else 
    Z_prime = (Z_r(i) - Z_r(i-1))/(r_r(i)-r_r(i-1));
    r = r - ((Z/Z_prime))
    r_r = [r_r;r];
end

end

% Saving steady state
V_st = V;
gg_st = gg;
K_st = K;
w_st = w;
r_st = r;
A_st = p.A;

% Transition dynamics 

gg0 = gg_st;

clear p.Delta r w p.A gg

%preallocation
gg = cell(N+1,1);
K_t = zeros(N,1);
K_out = zeros(N,1);
r_t = zeros(N,1);
w_t = zeros(N,1);

%initial guess
K_t = K_st*ones(N,1);

for it=1:p.maxit


w_t = (1-p.alpha)*A_t.*K_t.^p.alpha*L^(-p.alpha);
r_t = p.alpha*A_t.*K_t.^(p.alpha-1)*L^(1-p.alpha) - p.d;

V = V_st;


for n=1:N
    
    % 4-1. Compute the derivative of the value function 
    dVf = Df*V;
    dVb = Db*V;
    
    % 4-2. Boundary conditions
    dVb(1,:) = p.mu(w_t(n)*zz(1,:) + r_t(n).*aa(1,:)); % a>=a_min is enforced (borrowing constraint)
    dVf(end,:) = p.mu(w_t(n)*zz(end,:) + r_t(n).*aa(end,:)); % a<=a_max is enforced which helps stability of the algorithm
    
    I_concave = dVb > dVf; % indicator whether value function is concave (problems arise if this is not the case)
    
    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf);
    cb = p.inv_mu(dVb);
    
    % 4-4. Compute the optimal savings
    sf = w_t(n)*zz + r_t(n).*aa - cf;
    sb = w_t(n)*zz + r_t(n).*aa - cb;
    
    % 4-5. Upwind scheme
    If = sf>0;
    Ib = sb<0;
    I0 = 1-If-Ib;
    dV0 = p.mu(w_t(n)*zz + r_t(n).*aa); % If sf<=0<=sb, set s=0
    
    dV_upwind = If.*dVf + Ib.*dVb + I0.*dV0;
    
    c = p.inv_mu(dV_upwind);
    
    % 4-6. Update value function: 
    % Vj^(n+1) = [(rho + 1/Delta)*I - (Sj^n*Dj^n+A_switch)]^(-1)*[u(cj^n) + 1/Delta*Vj^n]
    
    V_stacked = V(:); % 2I*1 matrix
    c_stacked = c(:); % 2I*1 matrix
    
    % A = SD
    SD_u = spdiags(If(:,1).*sf(:,1), 0, p.I, p.I)*Df + spdiags(Ib(:,1).*sb(:,1), 0, p.I, p.I)*Db; % I*I matrix
    SD_e = spdiags(If(:,2).*sf(:,2), 0, p.I, p.I)*Df + spdiags(Ib(:,2).*sb(:,2), 0, p.I, p.I)*Db; % I*I matrix
    SD = [SD_u, sparse(p.I, p.I);
    sparse(p.I, p.I), SD_e]; % 2I*2I matrix
    
    % P = A + A_switch
    P = SD + A_switch;
    P_t{n} = P;
    
    % B = [(rho + 1/Delta)*I - P]
    B = (p.rho + 1/p.Delta)*speye(2*p.I) - P; 
    
    % b = u(c) + 1/Delta*V
    b = p.u(c_stacked) + (1/p.Delta)*V_stacked;
    
    % V = B\b;
    V_update = B\b; % 2I*1 matrix
    V_change = V_update - V_stacked;
    V = reshape(V_update, p.I, 2); % I*2 matrix
    % V = .9*V + .1*reshape(V_update, p.I, 2); % I*2 matrix
    
    
    % % 3-6. Convergence criterion
    % dist(n) = max(abs(V_change));
    % if dist(n)<p.tol
    %     disp('Value function converged. Iteration = ')
    %     disp(n)
    %     break
    % end
end

gg{1}=gg0;
for n=1:N
    PT=P_t{n}';
    %Implicit method in Updating Distribution.
    gg_pre= (speye(2*p.I) - PT*dt)\reshape(gg{n}, [],1);
    gg{n+1}= reshape(gg_pre, [], 2);
    
    %check(n) = gg(:,n)'*ones(2*p.I,1)*da;
    % K_out(n) = gg{n}(1:p.I)'*a*da + gg{n}(p.I+1:2*p.I)'*a*da;
    K_out(n) = gg{n}(1:p.I)*a*da + gg{n}(p.I+1:2*p.I)*a*da;
    
end

dist_it(it) = max(abs(K_out - K_t));

dist_it(it)


K_t = .5.*K_out +(1-.5).*K_t;

if dist_it(it)<p.tol
    disp('Equilibrium Found')
    break
end

end

% TFP Sequence

set(gca, 'FontSize', 18)
plot(time, A_t, 'LineWidth', 2)
grid
xlabel('Time, t','FontSize', 14)
ylabel('TFP, A_t','FontSize', 14)
xlim([min(time) max(time)])

% Wage transition path

set(gca, 'FontSize', 18)
plot(time, w_t, 'LineWidth', 2)
grid
xlabel('Time, t','FontSize', 14)
ylabel('Wage, w_t','FontSize', 14)
xlim([min(time) max(time)])

% Interest rate transition path

set(gca, 'FontSize', 18)
plot(time, r_t, 'LineWidth', 2)
grid
xlabel('Time, t','FontSize', 14)
ylabel('Interest rate, r_t','FontSize', 14)
xlim([min(time) max(time)])

% Capital transition path

set(gca, 'FontSize', 18)
plot(time, K_t, 'LineWidth', 2)
grid
xlabel('Time, t','FontSize', 14)
ylabel('Capital, K_t','FontSize', 14)
xlim([min(time) max(time)])




