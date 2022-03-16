function obj = singleShootingObjective(ubo, idxbo, x0, tspan, U, D, p, ...
    scalingFactor, objectiveFunction, simulationModel, outputModel, simulationMethod, opts)
% SINGLESHOOTINGOBJECTIVE evaluate the objective function in a dynamical
% optimization problem for determining the optimal insulin bolus.
%
% SYNOPSIS:
%   obj = singleShootingObjective(ubo, idxbo, x0, tspan, U, D, p, ...
%     objectiveFunction, simulationModel, outputModel)
%
% DESCRIPTION:
% This function evaluates the objective function in a dynamic optimization
% problem that is transcribed using the single-shooting approach, i.e.,
% given the initial condition and the manipulated inputs, the initial value
% problem, which constitutes the dynamical constraints, is solved for the
% states and the objective function is evaluated using these states.
% 
% Here, it is specifically assumed that the decision variables is an
% insulin bolus which is given during the first control interval. The
% objective function, the simulation model, the initial condition, the
% disturbances, and the model parameters are provided by the user.
%
% REQUIRED PARAMETERS:
%   ubo         - the insulin bolus                             (scalar)
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
%   obj         - the value of the objective function
%
% DEPENDENCIES:
%
% See also computeOptimalBolus
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

% Sizes of control intervals
Ts = diff(tspan);

% Verify that all control intervals are of the same size
allControlIntervalsAreEquallySized = (min(Ts) == max(Ts));

if(~allControlIntervalsAreEquallySized)
    error('The control intervals need to be the same size!');
end

% Meal and meal bolus
U(idxbo, 1) = ubo;

% Simulate
[T, X] = openLoopSimulation(x0, tspan, U, D, p, simulationModel, simulationMethod, opts);

% Evaluate the outputs
Z = outputModel(X, p);

% Evaluate the objective function
obj = scalingFactor*objectiveFunction(T, Z, p);