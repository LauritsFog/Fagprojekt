function [T, X] = odeEulersExplicitMethodFixedStepSize ...
    (fun, tspan, x0, opts, varargin) %#ok
% ODEEULERSEXPLICITMETHODFIXEDSTEPSIZE Use Euler's explicit method with
% fixed step size to approximate the solution to an initial value problem.
%
%   [T, X] = odeEulersExplicitMethodFixedStepSize(simModel, tspan, x0, opt)
%
% DESCRIPTION:
% Approximate the solution to the initial value problem
% 
%   x(t0) = x0,
%   dx(t)/dt = f(t, x(t))
% 
% using Euler's explicit method.
% 
% Note that the right-hand side function, f, can depend on more quantities
% than the time and states.
%
% REQUIRED PARAMETERS:
%   fun    - function handle for evaluating the right-hand side function
%   tspan  - points in time where the solution is approximated (dimension: N+1   )
%   x0     - initial states                                    (dimension: nx    )
%   opts   - options structure
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   T - boundaries of control intervals (=tspan)    (dimension: N+1)
%   X - the states in the simulation model          (dimension: N+1 x nx)
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

% Number of time steps
N = numel(tspan) - 1;

% Number of states
nx = numel(x0);

% Allocate memory
T = zeros(N+1,  1);
X = zeros(N+1, nx);

% Initial time
tk = tspan(1);

% Copy initial condition
xk = x0;

% Store initial condition
T(1   ) = tspan(1);
X(1, :) = x0;

for k = 1:N
    % Next time
    tkp1 = tspan(k+1);
    
    % Time step size
    dtk = tkp1 - tk;
    
    % Evaluate the right-hand side function
    fk = feval(fun, tk, xk, varargin{:});
    
    % Compute the states at the next time step
    xkp1 = xk + fk*dtk;
    
    % Store solution
    T(k+1   ) = tkp1;
    X(k+1, :) = xkp1;
end