function f = mvpModel(t, x, u, d, p) %#ok
% MVPMODEL Evaluate the right-hand side function of the Medtronic Virtual
% Patient (MPV) model.
%
% SYNOPSIS:
%   f = mvpModel(t, x, u, d, p)
%
% DESCRIPTION:
% This model describes the dynamics of the blood glucose concentration in
% response to meals and subcutaneous insulin.
%
% REQUIRED PARAMETERS:
%   t - time
%   x - a vector of state variables       (dimension:  7)
%   u - a vector a manipulated inputs     (dimension:  2)
%   d - a vector of disturbance variables (dimension:  1)
%   p - a vector parameters               (dimension: 10)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   f - the right-hand side function    (dimension:  7)
%
% DEPENDENCIES:
%
% See also mvpOutput
%          generateMVPParameters
%          mvpModelSteadyStateWrapper
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

% Conversion factors
g2mg = 1e3; % Convert from gram to milligram

% Meal subsystem
D1  = x(1); % [g CHO]
D2  = x(2); % [g CHO]

% Insulin subsystem
Isc = x(3); % [mU/L] Subcutaneous insulin concentration
Ip  = x(4); % [mU/L] Plasma insulin concentration

% Glucose subsystem
Ieff = x(5); % [1/min]  Insulin effect
G    = x(6); % [mg/dL] Blood glucose concentration

% Glucose concentration measurement
Gsc  = x(7); % [mg/dL] Subcutaneous glucose concentration

% Manipulated inputs
uba = u(1); % [ mU/min] Insulin basal flow rate
ubo = u(2); % [ mU/min] Insulin bolus flow rate

% Parameters
tau1    = p( 1); % [min]            Insulin absorption time
tau2    = p( 2); % [min]            Insulin absorption time
CI      = p( 3); % [dL/min]         Insulin clearance rate
p2      = p( 4); % [1/min]          Inverse of insulin action time constant
SI      = p( 5); % [(dL/mU)/min]    Insulin sensitivity
GEZI    = p( 6); % [1/min]          Glucose effectiveness
EGP0    = p( 7); % [(mg/dL)/min]    Endogenous glucose production
VG      = p( 8); % [dL]             Glucose distribution volume
taum    = p( 9); % [min]            Meal time constant
tausc   = p(10); % [min]            Subcutaneous insulin time constant

% Meal rate of appearance
RA = g2mg*D2/(VG*taum); % [(mg/dL)/min]

% Allocate memory
f = zeros(7, 1);

% Meal compartments
f(1) =  d  - D1 /taum;
f(2) = (D1 - D2)/taum;

% Insulin subsystem
f(3) = ((uba + ubo)/CI - Isc)/tau1;
f(4) = (Isc            - Ip )/tau2;

% Glucose subsystem
f(5) = p2*(SI*Ip - Ieff);
f(6) = -(GEZI + Ieff)*G + EGP0 + RA;

% Subcutaneous glucose concentration (measured by continuous glucose
% monitors - CGMs)
f(7) = (G - Gsc)/tausc;