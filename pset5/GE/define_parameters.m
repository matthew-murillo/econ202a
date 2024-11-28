function p = define_parameters()

% This function defines the parameters needed for the Huggett_PE.m script

%% Economic Parameters
    
    % Relative risk aversion
    p.gamma = 2;

    % Discount rate
    p.rho = 0.05;

    % Exogenous interest rate
    p.r = 0.035;

    % Exongeous wage
    p.w = 1;


    % Depreciation rate
    p.d = 0.05;

    % Capital share
    p.alpha = 1/3;

    % Exogenous TFP
    p.A = 0.1;

    % Idiosyncratic productivity
    p.z_u = 1;
    p.z_e = 2;
    p.z = [p.z_u, p.z_e];
    p.zz = [p.z_u, p.z_e];

    % Transition rates
    p.lambda_u = 1/3;
    p.lambda_e = 1/3;
    p.lambda = [p.lambda_u, p.lambda_e];
    
%% Economic Functions
    
    % Utility funtion
    p.u = @(c) c.^(1-p.gamma)/(1-p.gamma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.gamma);

    % FOC: mu(c)=dV -> c=inv_mu(dV)
    p.inv_mu = @(dV) dV.^(-1/p.gamma);

%% Grid Parmaters

    p.kmin = 0;
    p.kmax = 20;

    % The number of grid points
    p.I = 500;

%% Tuning parameters

    % Step size: can be arbitrarily large in implicit method
    p.Delta = 1000;

    % The maximum number of value function iterations
    p.maxit = 1000;

    % Tolerance for value function iterations
    p.tol = 10^(-6);

    % Max number of iterations for interest rate
    p.Nr = 1000;

    p.tol_S = 10^(-5);

end