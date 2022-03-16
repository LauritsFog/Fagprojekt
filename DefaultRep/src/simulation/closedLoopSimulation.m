function [T, X, Y, U, ctrlState] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState0, simMethod, opts)
% CLOSEDLOOPSIMULATION Simulate a closed-loop control algorithm.
%
% SYNOPSIS:
%   [T, X, U] = closedLoopSimulation(x0, tspan, D, p, ...
%     simModel, observationModel, ctrlAlgorithm, ...
%     ctrlPar, ctrlModel, ctrlModelPar, ctrlState)
%
% DESCRIPTION:
% Perform a closed-loop simulation of a model-based control algorithm for
% given initial condition, control intervals, disturbance variables,
% parameters, simulation model, observation model, and control algorithm
% (including hyperparamters, control model, control model paramters, and
% initial controller state).
%
% REQUIRED PARAMETERS:
%   x0                  - initial state                                     (dimension: nx    )
%   tspan               - boundaries of the control intervals               (dimension: N+1   )
%   D                   - disturbance variables for each control interval   (dimension: nd x N)
%   p                   - parameters                                        (dimension: np    )
%   simModel            - simulation model          (function handle)
%   observationModel    - observation model         (function handle)
%   ctrlAlgorithm       - control algorithm         (function handle)
%   ctrlPar             - controller parameters
%   ctrlState0          - initial controller state                          (dimension: nc)
%   simMethod           - simulation method         (function handle)
%   opts                - an options structure (must contain the field Nk)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   T - boundaries of control intervals (=tspan)    (dimension:      N+1)
%   X - the states in the simulation model          (dimension: nx x N+1)
%   X - the observed variables                      (dimension: ny x N+1)
%   U - the computed manipulated inputs             (dimension: nu x N  )
%   ctrlState - matrix of controller states         (dimension: nc x N+1)
%
% DEPENDENCIES:
%
% See also 
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

% Initial time
t0 = tspan(1);

% Observed variables
y0 = observationModel(t0, x0, p);

% Determine the number of manipulated inputs
uDummy = ctrlAlgorithm(t0, NaN, NaN, ctrlPar, ctrlState0);

% Number of states and manipulated inputs
nx = numel(x0);
ny = numel(y0);
nu = numel(uDummy);
nc = numel(ctrlState0);

% Number of control intervals
N = numel(tspan)-1;

% Number of time steps in each control interval
Nk = opts.Nk;

% Allocate memory
T = zeros( 1, N+1);
X = zeros(nx, N+1);
Y = zeros(ny, N+1);
U = zeros(nu, N  );
ctrlState = zeros(nc, N+1);

% Store initial condition
T(   1) = t0;
X(:, 1) = x0;
Y(:, 1) = y0;
ctrlState(:, 1) = ctrlState0;

% Copy initial condition
tk = t0;
xk = x0;
yk = y0;

for k = 1:N
    % Times
    tk    = tspan(k  );
    tkp1  = tspan(k+1);
    
    % Controller state
    ctrlStatek = ctrlState(:, k);
    
    % Disturbance variables
    dk = D(:, k);
    
    % Compute manipulated inputs
    [uk, ctrlStatekp1] = ctrlAlgorithm(tk, yk, dk, ctrlPar, ctrlStatek);
    
    % Time interval
    tspank = linspace(tk, tkp1, Nk+1);
    
    % Solve initial value problem
    [Tk, Xk] = simMethod(simModel, tspank, xk, [], uk, dk, p);
    
    % States at the next time step
    tkp1 = Tk(end   );
    xkp1 = Xk(end, :)';
    
    % Observed variables at the next time step
    ykp1 = observationModel(tkp1, xkp1, p);

    % Store solution
    T(   k+1) = tkp1;
    X(:, k+1) = xkp1;
    Y(:, k+1) = ykp1;
    U(:, k  ) = uk;
    ctrlState(:, k+1) = ctrlStatekp1;
    
    % Update initial condition
    tk = tkp1;
    xk = xkp1;
    yk = ykp1;
end