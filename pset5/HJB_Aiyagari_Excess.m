%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: Aiyagari Asset Supply and Demand
% 
% Author: Matt Murillo
%
%
% Reference: Andreas Schaab, Kiyea Jin, aiyagari_poisson_asset_supply.m by Benjamin Moll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();

S_r = [];

Ir = 100;
rmin = -0.049;
rmax = 0.049;
rgrid = linspace(rmin,rmax,Ir);

for r = rgrid
 
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
S_r = [S_r;S];
end

% Asset supply and asset demand
rrr = linspace(-0.06,0.06,Ir);
K = (p.alpha*p.A./(max(rrr+p.d,0))).^(1/(1-p.alpha))*L; % maybe problems here
K_grid = linspace(p.amin,max(S_r),Ir);

set(gca,'FontSize',14)
hold on;
plot(S_r, rgrid, 'LineWidth', 2, 'Color', 'r')
plot(K, rrr, 'LineWidth', 2,'Color', 'b')
plot(K_grid, repmat(p.rho - .001, Ir), '--', 'Color', 'black')
plot(K_grid, repmat(-p.d + .001, Ir), '--', 'Color', 'black')
hold on;
ylabel('Interest rate, r','FontSize',16)
xlabel('Capital, K','FontSize',16)
ylim([min(rgrid)-.01 p.rho+.01])
maxS = max(S_r);
xlim([p.amin-0.01 maxS])
hold off

% Excess demand
SD = K' - S_r;
set(gca,'FontSize',14)
plot(rgrid, SD, 'LineWidth', 2, 'Color', 'r')
ylabel('Excess demand')
xlabel('Interest rate, r')
ylim([min(SD), 10])
xlim([min(rgrid), max(rgrid)])


