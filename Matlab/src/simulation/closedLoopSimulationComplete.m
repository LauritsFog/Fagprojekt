function [T, X, Y, U, ctrlState] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState0, simMethod, tzero, haltingiter, idxbo, rampingfunction, ...
    dg, dt, gridTime, opts)
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

% Integration halt time during bolus titration
tpause = 0;

% Determine the number of manipulated inputs
uDummy = ctrlAlgorithm(t0, NaN, NaN, ctrlPar, ctrlState0, tzero, tpause, haltingiter, rampingfunction);

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

% Initialize GRID
GRID = zeros(1,gridTime);

Dest = zeros(1,length(tspan));

% Insuling to carbohydrate ratio
ICR = 1/2;

% If 1 meal size estimation is used, otherwise grid algo is used
mealMethod = 0;

% Minimim meal sizes allowed
minMeal = 0.3;

% bolus scalar
alpha = 0.2;
beta = 0.5;
gamma = 0.1;

% Parameters for lowpassfilter
dtlowpass = 3;
tau = 6;

for k = 1:N
    % Times
    tkp1 = tspan(k+1);
    
    % Controller state
    ctrlStatek = ctrlState(:, k);
    
    % Disturbance variables
    dk = D(:, k);
    
    % Detecting meals
    if mealMethod == 0 % Using grid algo
        if k > gridTime
            [~,~,GRID]=GridAlgo(Y(k-gridTime:k),dg,dt,[],tspan(k-gridTime:k));
        end

        % If meal is detected and no meal has been detected for in past
        % gridTime
        if k > gridTime && nnz(GRID) > 0 && nnz(Dest(k-gridTime:k)) == 0
            dkest = 1;
        else
            dkest = 0;
        end

        Dest(k) = dkest;

        % If meal is detected, halt integration for some iterations and compute
        % optimal bolus
        if dkest ~= 0
            Ylowpass = lowpassfilter(Y(k-gridTime:k),dtlowpass,tau,gridTime);

            dYlowpass = Ylowpass(2:end)-Ylowpass(1:end-1);

            dYmean = mean(dYlowpass);

            tpause = haltingiter;
            
            ubok = alpha*Ylowpass(end)*haltingiter*max([0,beta+gamma*dYmean]);
        else
            ubok = 0;
        end
        
        if yk < ctrlPar(5)
            ubok = 0;
        end
    else % Using meal estimation
        if k > gridTime
            [~,MealEst,~,~] = MealSize(Y(k-gridTime:k),tspan(k-gridTime:k));
            
            if sum(MealEst) > minMeal
                dkest = sum(MealEst);
            else
                dkest = 0;
            end
        else
            dkest = 0;
        end
        
        Dest(k) = dkest;
        
        if dkest ~= 0 && tpause == 0
            tpause = haltingiter;
        end
        
        ubok = ICR*dkest;
    end
    
    % Decrement tpause until 0
    if tpause ~= 0
        tpause = tpause - 1;
    end
    
    % Compute manipulated inputs
    [uk, ctrlStatekp1] = ctrlAlgorithm(tk, yk, dk, ctrlPar, ctrlStatek, tzero, tpause, haltingiter, rampingfunction);
    
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