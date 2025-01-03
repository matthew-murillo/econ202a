function p = define_parameters()

% This function defines the parameters needed for the HJB_ramsey_implicit.m script

%% Economic Parameters

    % Discount rate
    p.rho = 0.05;

    % Relative risk aversion coefficient
    p.sigma = 2;
    
    % Interest rate;
    p.r = 0.03

    % Transition rates
    p.lambda = [1/3, 1/3]

    % Income
    p.z = [0.2, 0.1]

    % Asset min
    p.amin = -0.02

    % Asset max (arbitrary)
    p.amax = 10

%% Economic Functions
    
    % Utility function
    p.u = @(c) (c.^(1-p.sigma))./(1-p.sigma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.sigma);

    % Inverse of marginal utility
    % Note: FOC: mu(c) = dv(a) -> c = inv_mu(dv)
    p.inv_mu = @(dv) dv.^(-1/p.sigma);

    % Production function
    p.f = @(k) p.A * k.^p.alpha;

     
      

%% Grid Paramters

    % The number of grid points
    p.I = 500;

    % Higher klim implies broader coverage
    p.klim = 1.5;

%% Tuning Parameters
    
    % The maximum number of iteration
    p.maxit = 1000;

    % Convergence criterion, and
    p.tol = 1e-8;

    % Step size (Delta)
    p.Delta = 1000;

end