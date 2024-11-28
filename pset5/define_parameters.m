function p = define_parameters()

% This function defines the parameters needed for the HJB_Huggett.m script

%% Economic Parameters

    % Discount rate
    p.rho = 0.05;

    % Relative risk aversion coefficient
    p.sigma = 2;


    % Income
    p.z_u = 1;
    p.z_e = 2;
    p.zz = [p.z_u, p.z_e];

    % Transition rates
    p.lambda_u = 1/3;
    p.lambda_e = 1/3;
    p.lambda = [p.lambda_u, p.lambda_e];

    % Capital share
    p.alpha = 1/3;

    % Productivity
    p.A = 0.1;

    % Depreciation
    p.d = 0.05;


    
%% Economic Functions
    
    % Utility function
        % if sigma == 1
        % p.u = @(c) log(c);
        % else
        % p.u = @(c) (c.^(1-sigma))./(1-sigma);
        % end
    p.u = @(c) (c.^(1-p.sigma))./(1-p.sigma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.sigma);

    % Inverse of marginal utility
    % Note: FOC: mu(c) = dv(a) -> c = inv_mu(dv)
    p.inv_mu = @(dv) dv.^(-1/p.sigma);

%% Grid Paramters

    % The lower bound of the state space (borrowing constraint)
    p.amin = 0;

    % The upper bound of the state space
    p.amax = 20;

    % The number of grid points
    p.I = 500;

%% Tuning Parameters
    
    % The maximum number of iteration
    p.maxit = 1000;

    % Convergence criterion
    p.tol = 1e-8;

    % Step size (Delta)
    p.Delta = 1000;

end