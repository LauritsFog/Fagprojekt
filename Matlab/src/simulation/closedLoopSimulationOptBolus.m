function [T, X, Y, U, ctrlState] = closedLoopSimulationOptBolus(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState0, simMethod, tzero, haltingiter, ...
    ubo0, idxbo, scalingFactor, objectiveFunction, outputModel, rampingfunction, opts)
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
% Svensson, J., Ranjan, A. G., N??rgaard, K. J??rgensen, J. B., 2022.
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
% Asbj??rn Thode Reenberg
% John Bagterp J??rgensen

% Initial time
t0 = tspan(1);

% Observed variables
y0 = observationModel(t0, x0, p);

% Integration halt time during bolus titration
tpause = 0;

% Time span used when computing optimal bolus
% tspanbolus = tspan(1:haltingiter);
tspanbolus = tspan;

% Determine the number of manipulated inputs
uDummy = ctrlAlgorithm(t0, NaN, NaN, ctrlPar, ctrlState0, tzero, tpause);

% Number of states and manipulated inputs
nx = numel(x0);
ny = numel(y0);
nu = numel(uDummy);
nc = numel(ctrlState0);

% Number of control intervals
N = numel(tspan)-1;

% Basal during bolus titration
% Ubolus = repmat([ctrlPar(6);0], 1, length(tspanbolus));

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
    tkp1 = tspan(k+1);
    
    % Controller state
    ctrlStatek = ctrlState(:, k);
    
    % Disturbance variables
    dk = D(:, k);
    
    % If meal is detected, halt integration for some iterations and compute
    % optimal bolus
    if dk ~= 0
        tpause = haltingiter;
        
        Ubolus = simulatePID(tk, xk, yk, dk, Nk, p, ctrlPar, ctrlStatek, ctrlAlgorithm, simModel, simMethod, observationModel, N, tzero, haltingiter, rampingfunction);
        
        % Ubolus = repmat([ctrlPar(6);0], 1, length(D));
        
        % plot(linspace(1,length(Ubolus),length(Ubolus)),Ubolus)
        
        Dtemp = zeros(1,length(D));
        Dtemp(1) = D(k);
        
        [ubok, flag] = computeOptimalBolus(ubo0, idxbo, xk, tspanbolus, Ubolus, Dtemp, p, ...
        scalingFactor, objectiveFunction, simModel, outputModel, simMethod, opts);
        
        % If fsolve did not converge, throw an error
        if(flag ~= 1), error ('fmincon did not converge!'); end
    else
        ubok = 0;
    end
    
    % Compute manipulated inputs
    [uk, ctrlStatekp1] = ctrlAlgorithm(tk, yk, dk, ctrlPar, ctrlStatek, tzero, tpause, haltingiter, rampingfunction);
    
    % Decrement tpause until 0
    if tpause ~= 0
        tpause = tpause - 1;
    end
    
    % Set optimal bolus 
    uk(idxbo) = ubok;
    
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