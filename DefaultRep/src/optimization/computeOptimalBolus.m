function [ubo, flag] = computeOptimalBolus(ubo0, idxbo, x0, tspan, U, D, p, ...
    scalingFactor, objectiveFunction, simulationModel, outputModel, simulationMethod, opts)
% COMPUTEOPTIMALBOLUS compute an optimal insulin bolus for a given
% objective function and simulation model.
%
% SYNOPSIS:
%   [ubo, Tbo, Xbo, flag] = computeOptimalBolus(ubo0, idxbo, x0, tspan, U, D, p, ...
%     objectiveFunction, simulationModel, outputModel)
%
% DESCRIPTION:
% This function uses fmincon to solve a dynamic optimization problem where
% the decision variable is the insulin bolus administered during the first
% control interval. The decision variable is chosen such as to minimize a
% user-provided objective function over a user-defined prediction horizon.
% The objective function, the simulation model, the output model, the
% initial condition, the disturbances, and the model parameters are
% provided by the user.
%
% REQUIRED PARAMETERS:
%   ubo0        - initial guess of the optimal insulin bolus    (scalar)
%   idxbo       - index of the insulin bolus                    (scalar)
%   x0          - initial condition                             (dimension: nx)
%   tspan       - vector of boundaries of the control intervals (dimension: N+1)
%   U           - matrix of manipulated inputs                  (dimension: nu x N)
%   D           - matrix of disturbances                        (dimension: nd x N)
%   p           - vector of parameters                          (dimension: np)
%   scalingFactor       - a scaling factor for the objective function (scalar)
%   objectiveFunction   - a function handle to the objective function
%   simulationModel     - a function handle to the simulation model
%   outputModel         - a function handle to the output model
%   simulationMethod    - a function handle to the simulation method/function
%   opts                - an options structure (must contain the field Nk)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   ubo         - the (locally) optimal insulin bolus
%   flag        - exit flag returned by fmincon
%
% DEPENDENCIES:
%
% See also singleShootingObjective
% 
% REFERENCES
% [1] Kanderian, S. S., Weinzimer, S., Voskanyan, G., and Steil, G. M.
% (2009). Identification of intraday metabolic profiles during closed-loop
% glucose control in individuals with type 1 diabetes. Journal of Diabetes
% Science and Technology, 3(5), 1047-1057.
% [2] Facchinetti, A., Favero, S. D., Sparacino, G., Castle, J. R., Ward,
% W. K., and Cobelli, C. (2014). Modeling the glucose sensor error. IEEE
% Transactions on Biomedical Engineering, 61(3), 620-629.
% [3] Reenberg, A. T., Ritschel, T. K. S., Lindkvist, E. B., Laugesen, C.,
% Svensson, J., Ranjan, A. G., Nørgaard, K. Jørgensen, J. B., 2022.
% Nonlinear Model Predictive Control and System Identification for a Dual-
% hormone Artificial Pancreas. In submission. arXiv: 2202.13938.
% 
% CONTACT INFORMATION
%  info@diamatica.com
%  tobk@dtu.dk
% athre@dtu.dk
%  jbjo@dtu.dk
% 
% AUTHORS
% Tobias K. S. Ritschel
% Asbjørn Thode Reenberg
% John Bagterp Jørgensen

% There are no equality constraints
A = [];
b = [];

% There are no inequality constraints (except bounds)
Aeq = [];
beq = [];

% Bounds
lb = 0;
ub = 1e9;

% There are no nonlinear constraints
nonlcon = [];

% Options
fminconOpts = optimoptions('fmincon', ...
    'Display', 'none');

% Solve the dynamic optimization problem
[ubo, ~, flag] = fmincon(@singleShootingObjective,  ...
    ubo0, A, b, Aeq, beq, lb, ub, nonlcon, fminconOpts,    ...
    idxbo, x0, tspan, U, D, p, scalingFactor,       ...
    objectiveFunction, simulationModel, outputModel, simulationMethod, opts);