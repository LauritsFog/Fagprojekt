function [uk, ctrlState] = pidController ...
    (tk, yk, dhatk, ctrlPar, ctrlState) %#ok
% PIDCONTROLLER PID controller for controlling the insulin flow rate.
%
% SYNOPSIS:
%   [uk, ctrlState] = pidController ...
%     (tk, yk, dhatk, ctrlPar, ctrlModel, ctrlModelPar, ctrlState)
%
% DESCRIPTION:
% This function implements a discretized proportional-integral-derivative
% (PID) controller for controlling the blood glucose concentration.
%
% REQUIRED PARAMETERS:
%   tk              - time
%   yk              - observed variables at time tk
%   dhatk           - estimated disturbance variables at time tk
%   ctrlPar         - controller parameters
%   ctrlModel       - control model
%   ctrlModelPar    - control model parameters
%   ctrlState       - controller state
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   uk          - a vector of manipulated inputs
%   ctrlState   - the updated controller state
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

% Unpack control parameters
Ts      = ctrlPar(1); % [min]    Sampling time
KP      = ctrlPar(2); %          Proportional gain
KI      = ctrlPar(3); %          Integrator gain
KD      = ctrlPar(4); %          Derivative gain
ybar    = ctrlPar(5); % [mg/dL]  Target blood glucose concentration
ubar    = ctrlPar(6); % [mU/min] Nominal insulin flow rate

% Unpack control state
Ikm1 = ctrlState(1); %           Value of integral at previous time step
ykm1 = ctrlState(2); % [mg/dL]   Previous observed glucose concentration

% Setpoint error
ek = yk - ybar;

% Derivative
dek = (yk - ykm1)/Ts;

% Proportional term
Pk = KP*ek;

% Integral term
Ik = Ikm1 + KI*ek*Ts;

% Derivative term
Dk = KD*dek;

% Basal insulin flow rate
ubak = ubar + Pk + Ik + Dk;

% Bolus insulin flow rate
ubok = 0;

% Manipulated inputs (must be non-negative)
uk = [ubak; ubok];

% Controller state
ctrlState = [Ik; yk];