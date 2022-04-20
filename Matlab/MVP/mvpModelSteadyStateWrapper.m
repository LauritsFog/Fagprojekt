function R = mvpModelSteadyStateWrapper(w, t, d, p, Gs)
% MVPMODELSTEADYSTATEWRAPPER Wrapper used to identify the steady state
% of the Medtronic Virtual Patient (MVP) model.
%
% SYNOPSIS:
%   R = mvpModelSteadyStateWrapper(w, t, d, p, Gs)
%
% DESCRIPTION:
% This function is used together with a root-finding algorithm, e.g., one
% implemented in Matlabs fsolve, in order to find the steady state of the
% Medtronic Virtual Patient (MVP) model for a given glucose concentration
% and zero insulin and glucagon boli:
% 
%   [f(t, w, d, p); A w - b]
%
% REQUIRED PARAMETERS:
%   w  - a vector of states and manipulated inputs  (dimension:  9)
%   t  - time
%   d  - a vector of disturbance variables          (dimension:  1)
%   p  - a vector parameters                        (dimension: 10)
%   Gs - the steady state blood glucose concentration [mg/dL]
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   R - the residual equations (dimension: 9)
%
% DEPENDENCIES:
% mvpModel
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

% Extract states
x = w(1:end-2);

% Extract manipulated inputs
u = w(end-1:end);

% Specified state variable
G = x(6); % [mg/dL]

% Specified manipulated input
ubo = u(2); % [mU/min]

% Evaluate the right-hand side of the Hovorka model
R1 = mvpModel(t, x, u, d, p);

% Specification equation for the blood glucose concentration
R2 = G - Gs; % [mg/dL]

% Specification equation for bolus insulin flow rate
R3 = ubo;

% Collect residual equations
R = [R1; R2; R3];