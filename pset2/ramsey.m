%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script: ramsey.m
%
% Author: Matt Murillo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

t = linspace(p.tmin, p.tmax, p.I)';
dt = (p.tmax-p.tmin)/(p.I-1);

%% 3. STEADY STATES Analytical
% kss = ((rho + theta*g)/(alpha*A))^(1/(alpha-1))
% css = kss * (rho + (theta*g) - (alpha*(n + g)))/alpha

kss = ((p.rho + p.theta * p.g) / (p.alpha * p.A))^(1 / (p.alpha - 1));
css = kss * (p.rho + (p.theta*p.g) - (p.alpha*(p.n + p.g)))/p.alpha;
rss = p.f_prime(kss);
%% 4. STEADY STATE FSOLVE


diff = @(x) [p.f(x(1)) - x(2) - (p.n + p.g) * x(1); p.f_prime(x(1)) - p.rho - p.theta * p.g]; 

k0 = 10;
c0 = 10;

x0 = [k0; c0];  % Initial guess for k and c

options = optimoptions('fsolve', 'Display', 'off');
x = fsolve(diff, x0, options);

% Steady state k and c and interest rate
kss_N = x(1);
css_N = x(2);
rss_N = p.f_prime(kss_N);
%% 5. SHOOTING ALGORITHM

% 4-1. Find the initial value of consumption that satisfies terminal boundary condition

tic;
% Objective function that calculates the difference between terminal capital k(T) and steady-state kss
% Note: k(T)-kss is a function of c0
diff = @(c0) terminal_condition(c0, p.k0, kss, p.f, p.f_prime, p.rho, p.theta, p.g, p.n, dt, p.I);

% Guess an initial value of consumption
c0_guess = 1;

% Use fsolve to find the initial consumption c0 that makes k(T) = k_ss
% Note: X = fsolve(FUN, X0, options) starts at the matrix X0 and tries to solve the equations in FUN.
% Set OPTIONS = optimoptions('fsolve','Algorithm','trust-region'), and then pass OPTIONS to fsolve.
options = optimoptions('fsolve', 'TolFun', p.tol, 'Display', 'iter');
c0 = fsolve(diff, c0_guess, options);

% 4-2. Forward simulate with the updated initial consumption

[k, c] = forward_simulate(c0, p.k0, p.f, p.f_prime, p.rho, p.theta, p.g, p.n, dt, p.I);
toc;

% 5-1. Evolution of capital and consumption and interest rate
rt = p.f_prime(k);

figure;
subplot(3,1,1);
plot(t, k, 'r-', 'LineWidth', 2);
xlabel('Time'); ylabel('Capital k(t)');
title('Capital Accumulation over Time');

subplot(3,1,2);
plot(t, c, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('Consumption c(t)');
title('Consumption Growth over Time');

subplot(3,1,3);
plot(t, rt, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('Interest rate r(t)');
title('Interest rate over Time');

saveas(gcf, "fig1.png")

%% 6 Newton's Method 

% Initial guess
k0 = 10;
c0 = 10;

b0 = [k0;c0]

for i = 1: 1000
kt= b0(1)
ct= b0(2)

J_inv = inv([(p.alpha * p.A * kt^(p.alpha-1)) - (p.n + p.g), -1; (ct*p.alpha*(p.alpha-1)*p.A*kt^(p.alpha-2))/p.theta, (p.alpha*p.A*kt^(p.alpha-1)-p.rho-(p.theta*p.g))/p.theta])

F = [p.A*kt^(p.alpha)-ct-((p.n + p.g)*kt); ct*((p.alpha*p.A*kt^(p.alpha-1)-p.rho-(p.theta*p.g))/p.theta)]

b1 = b0 - J_inv*F

if norm(b1-b0) < p.tol
        break;
end

b0 = b1

end


kss_J = b1(1);
css_J = b1(2);
rss_J = p.f_prime(kss_J);


