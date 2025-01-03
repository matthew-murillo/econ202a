%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_ramsey_implicit_upwind
% 
% Author: Matt Murillo
%
% Description:
% This MATLAB script implements implicit method to solve the HJB equation
% of the deterministic Ramsey Growth Model using upwind scheme.
%
% Reference:
% HJB_NGM_implicit.m by Benjamin Moll
% ramsey_implicit.m by Pontus Rendahl
% HJB_ramsey_implicit_upwind.m by Kiyea Jin

% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. PRE-ITERATION INITIALIZATION
% 4. VALUE FUNCTION ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters3();

%% 2. INITIALIZE GRID POINTS
% Grid of assets
a = linspace(p.amin, p.amax, p.I)';
% Adjusting for vector 

% Initial guess for each a
v0 = p.u(p.z + p.r*a)/p.rho;
V = v0;

% Stacking the vector
V = reshape(v0, [], 1);

da = (p.amax - p.amin)/(p.I-1);

%% 3. PRE-ITERATION INITIALIZATION

% 3-1. Construct the forward and backward differential operator 
% Df such that Df*V=dVf and Db such that Db*V=dVb. 

% This matrix will be used to get (V(t+1) - V(t))/da (the finite forward difference
% approximiation of the derivative of the value function)

Df = zeros(p.I, p.I);
for i = 1:p.I-1
    Df(i,i) = -1/da; Df(i,i+1) = 1/da;
end
% Adjust to accomdate long V
Df = blkdiag(Df,Df);
% This matrix will be used to get (V(t) - V(t-1))/da (the finite backward difference
% approximiation of the derivative of the value function)
Db = zeros(p.I, p.I);
for i = 2:p.I
    Db(i,i-1) = -1/da; Db(i,i) = 1/da;
end
Db = blkdiag(Db,Db);

% Constructing the transition matrix

% Transition rates
A1 = speye(p.I).*-p.lambda(1); % Rate of leaving employment
A2 = speye(p.I).*p.lambda(1); % Rate of becoming employed

A3 = speye(p.I).*p.lambda(2); % Rate of becoming unemployed
A4 = speye(p.I).*-p.lambda(2); % Rate of leaving unemployment

A = [A1,A2;A3,A4];

% 3-3. Pre-allocate arrays for solutions

dVf = zeros(p.I,1);
dVb = zeros(p.I,1);
c = zeros(p.I,1);

%% 4. VALUE FUNCTION ITERATION

tic;

for n = 1:p.maxit

    % 4-1. Compute the derivative of the value function
        dVf = Df*V;
        dVb = Db*V;

    % BOUNDARY CONDITIONS
        boundary_f = p.mu(p.z + p.r.*p.amax); % Forward boundary condtion
        boundary_b = p.mu(p.z + p.r.*p.amin); % Backward boundary condition

        dVf(dVf == 0) = boundary_f;
        dVb(dVb == 0) = boundary_b;
        % dVf(end) = p.mu(p.f(k_max) - p.delta*k_max); % k<=k_max is enforced which helps stability of the algorithm
        % dVb(1)   = p.mu(p.f(k_min) - p.delta*k_min); % k>=k_min is enforced which helps stability of the algorithm

    % 4-2. Compute the optimal consumption
        cf = p.inv_mu(dVf);
        cb = p.inv_mu(dVb);

    % 4-3. Compute the optimal savings
        sf = reshape(repmat(p.z, p.I, 1), [],1) + repmat(p.r*a, 2, 1) - cf;
        sb = reshape(repmat(p.z, p.I, 1), [],1) + repmat(p.r*a, 2, 1) - cb;
      
   
    % UPWIND SCHEME
        If = sf>0;
        Ib = sb<0;
        I0 = 1-If-Ib;

        dV0 = p.mu(reshape(repmat(p.z, p.I, 1), [],1) + repmat(p.r*a, 2, 1));
        dV_upwind = dVf.*If + dVb.*Ib + dV0.*I0;

        c = p.inv_mu(dV_upwind);

    % 4-4. Update the value function: V^(n+1) = [(rho+1/Delta)*I - SD]^(-1)[u(c) + 1/Delta*V^n]
    
        % B = [(rho+1/Delta)*I - SD]
        Sf = diag(sf.*If);
        Sb = diag(sb.*Ib);
        SD = Sf*Df + Sb*Db;
        B_mat = (p.rho + 1/p.Delta)*eye(p.I);
        B = blkdiag(B_mat, B_mat) - (SD + A);

        % b = [u(c) + 1/Delta*V^n]
        b = p.u(c) + 1/p.Delta*V;

        % V^(n+1) = B^(-1)*b
        V_update = B\b;
        
        % Update the value function
        V_change = V_update - V;
        V = V_update;

    % 4-5. Check convergence
          
        dist(n) = max(abs(V_change));

        if dist(n)<p.tol
        disp('Value function converged. Iteration = ')
        disp(n)
        break
        end

end

% Optimal saving
s = reshape(repmat(p.z, p.I, 1), [],1) + repmat(p.r*a, 2, 1) - c;


c = reshape(c, [], 2);
s = reshape(s, [], 2);
V = reshape(V, [], 2);
% 
% set(gca,'FontSize',14)
% plot(a,c,'LineWidth',2)
% grid
% xlabel('a')
% ylabel('c_i(a)')
% xlim([p.amin p.amax])
% legend('Employed','Unemployed', 'Location', 'northwest')
% saveas(gcf, 'fig1.png')
% 
% set(gca,'FontSize',14)
% plot(a,s,'LineWidth',2)
% grid
% xlabel('a')
% ylabel('s')
% xlim([p.amin p.amax])
% legend('Employed','Unemployed', 'Location', 'northwest')
% saveas(gcf, 'fig2.png')
% 
% set(gca,'FontSize',14)
% plot(a,V,'LineWidth',2)
% grid
% xlabel('a')
% ylabel('V')
% xlim([p.amin p.amax])
% legend('Employed','Unemployed', 'Location', 'northwest')
% saveas(gcf, 'fig3.png')
% 
