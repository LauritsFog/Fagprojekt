function [T, X] = openLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, opts)
% OPENLOOPSIMULATION Open-loop simulation
%
%   [T, X] = openLoopSimulation(x0, tspan, U, D, p, simModel, simMethod)
%
% DESCRIPTION:
% Perform an open-loop simulation for given initial condition, control
% intervals, disturbance variables, parameters, and simulation model.
%
% REQUIRED PARAMETERS:
%   x0                  - initial state                                     (dimension: nx    )
%   tspan               - boundaries of the control intervals               (dimension: N+1   )
%   U                   - disturbance variables for each control interval   (dimension: nu x N)
%   D                   - disturbance variables for each control interval   (dimension: nd x N)
%   p                   - parameters                                        (dimension: np    )
%   simModel            - simulation model          (function handle)
%   simMethod           - simulation method         (function handle)
%   opts                - an options structure (must contain the field Nk)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   T - boundaries of control intervals (=tspan)    (dimension:      N+1)
%   X - the states in the simulation model          (dimension: nx x N+1)
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

% Number of control intervals
N = numel(tspan)-1;

% Number of time steps in each control interval
Nk = opts.Nk;

% Number of states
nx = numel(x0);

% Allocate memory
T = zeros( 1, N+1);
X = zeros(nx, N+1);

% Initial condition in each control interval
xk = x0;

% Store solution
T(   1) = tspan(1);
X(:, 1) = x0;

for k = 1:N
    % Times
    tk    = tspan(k  );
    tkp1  = tspan(k+1);
    
    % Manipulated inputs and disturbance variables
    uk = U(:, k);
    dk = D(:, k);
    
    % Time interval
    tspank = linspace(tk, tkp1, Nk+1);
    
    % Solve initial value problem
    [Tk, Xk] = simMethod(simModel, tspank, xk, [], uk, dk, p);
    
    % Update initial condition
    xk = Xk(end, :)';
    
    % Store solution
    T(   k+1) = Tk(end   )';
    X(:, k+1) = Xk(end, :)';
end