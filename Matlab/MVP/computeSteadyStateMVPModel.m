function [xs, us, flag] = computeSteadyStateMVPModel(t, p, Gs)
% COMPUTESTEADYSTATEMVPMODEL Compute the steady state of the Medtronic
% Virtual Patient (MVP) model.
%
% SYNOPSIS:
%   [xs, us, flag] = computeSteadyStateMVPModel(t, p, Gs)
%
% DESCRIPTION:
% Solve a combination of the nonlinear right-hand side function and a set
% of linear specification equations for the steady state of the Medtronic
% Virtual Patient (MVP) model at given blood glucose concentration and zero
% insulin and glucagon boli:
% 
%   f(t, w, d, p) = 0,
%   A(p) w = b(p).
% 
% Here, w = [x; u]. This set of nonlinear and linear equations are solved
% using fsolve to a hardcoded tolerance of 10^(-10).
% 
% Note that for the MVP model, the steady state can also be computed
% analytically. However, that is not the case for more complicated models.
%
% REQUIRED PARAMETERS:
%   t  - time
%   p  - a vector parameters               (dimension: 10)
%   Gs - the steady state blood glucose concentration [mg/dL]
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   xs   - the steady state                    (dimension: 7)
%   us   - the steady state manipulated inputs (dimension: 2)
%   flag - the flag (integer) returned by fsolve
%
% DEPENDENCIES:
% mvpModel
% mvpModelSteadyStateWrapper
%
% See also mvpModel
%          mvpOutput
%          generateMVPParameters
%          computeSteadyStateMVPModel
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

% fsolve options
opts = optimoptions('fsolve',       ...
    'OptimalityTolerance',  1e-10,  ...
    'Display',              'off');

% The disturbance variables
d = 0;

% Initial guess of the states
xs0 = naturalSteadyStateMVPModel(p);

% Initial guess of the manipulated inputs
us0 = [0; 0];

% Initial guess
ws0 = [xs0; us0];

% Solve for the steady state
[ws, ~, flag] = fsolve(@mvpModelSteadyStateWrapper, ...
    ws0, opts, t, d, p, Gs);

% Steady state
xs = ws(1:end-2);

% Steady state manipulated inputs
us = ws(end-1:end);