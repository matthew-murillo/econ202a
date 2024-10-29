%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script: euler.m
%
% Author: Matt Murillo

% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. ANALYTICAL SOLUTION: INTEGRATING FACTOR METHOD
% 4. NUMERICAL SOLUTION: FINITE DIFFERENCE METHOD
% 5. PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;
clc;

%% 1. DEFINE PARAMETERS

C_T = 2;
r_0 = 0.05
alpha = 0.01
theta = 2
rho = 0.03
T = 10
T_steps = 100

%% 2. INITIALIZE GRID POINTS

t = linspace(0, T, T_steps)';
dt = (max(t)-min(t))/(T_steps);


%% 3. ANALYTICAL SOLUTION: INTEGRATING FACTOR METHOD

% 3-1. Pre-allocate arrays for solutions
C_analytical = zeros(T_steps,1);

% 3-2. Analytical solution
for i=1:T_steps
   
    C_analytical(i) = C_T * exp(-(((rho-r_0)/theta)*t(i) - (alpha/(2*theta))*t(i)^2)) * exp(((rho-r_0)/theta)*T - ((alpha/(2*theta))*T^2))
    
end


%% 4. NUMERICAL SOLUTION: FINITE DIFFERENCE METHOD

% 4-1. Pre-allocate arrays for solutions
C_numerical  = zeros(T_steps,1);
C_numerical(end) = C_T;

% 4-2. Numerical solution

for i=1:1:T_steps-1
    
    C_numerical(end - i) = C_numerical(end-(i-1)) * (1 + dt*(r_0 + alpha*t(end- (i-1)) - rho)/(theta))^-1
    
end

C_numerical_100 = C_numerical;

T_steps = 10
t_10 = linspace(0, T, T_steps)';
dt = (max(t_10)-min(t_10))/(T_steps);

% 4-1. Pre-allocate arrays for solutions
C_numerical  = zeros(T_steps,1);
C_numerical(end) = C_T;

% 4-2. Numerical solution

for i=1:1:T_steps-1
    
    C_numerical(end - i) = C_numerical(end-(i-1)) * (1 + dt*(r_0 + alpha*t_10(end- (i-1)) - rho)/(theta))^-1
    
end

C_numerical_10 = C_numerical;


%% 5. PLOT
% Plot both solutions for comparison
figure;
plot(t, C_analytical, 'r-', 'LineWidth', 2); hold on;
xlabel('Time', 'FontSize', 11);
ylabel('Consumption (C)', 'FontSize', 11);
legend('Analytical Solution','FontSize', 12);
title('Analytical Solution', 'FontSize', 15);
grid on;
saveas(gcf, "fig1.png")
clf


% Plot both solutions for comparison
figure;
plot(t, C_analytical, 'r-', 'LineWidth', 2); hold on;
plot(t, C_numerical_100, 'bo--', 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 11);
ylabel('Consumption (C)', 'FontSize', 11);
legend('Analytical Solution', 'Numerical Solution', 'FontSize', 12);
title('Comparison of Analytical and Numerical Solutions (Time step: 100)', 'FontSize', 15);
grid on;
saveas(gcf, "fig2.png")
clf


% Plot both solutions for comparison
figure;
plot(t, C_analytical, 'r-', 'LineWidth', 2); hold on;
plot(t_10, C_numerical_10, 'bo--', 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 11);
ylabel('Consumption (C)', 'FontSize', 11);
legend('Analytical Solution', 'Numerical Solution', 'FontSize', 12);
title('Comparison of Analytical and Numerical Solutions (Time step: 10)', 'FontSize', 15);
grid on;
saveas(gcf, "fig3.png")
clf
