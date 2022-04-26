function [xs, us] = naturalSteadyStateMVPModel(p)
% NATURALSTEADYSTATEMVPMODEL Compute the natural steady state of the
% Medtronic Virtual Patient (MVP) model.
%
% SYNOPSIS:
%   [xs, us] = naturalSteadyStateMVPModel(p)
%
% DESCRIPTION:
% Compute the steady state of the Medtronic Virtual Patient (MVP) model
% corresponding to no meals or insulin or glucagon derived analytically.
%
% REQUIRED PARAMETERS:
%   p - a vector parameters (dimension: 31)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   xs   - the steady state                    (dimension: 16)
%   us   - the steady state manipulated inputs (dimension:  3)
%
% DEPENDENCIES:
%
% See also generateExtendedHovorkaParameters
%          extendedHovorkaModelSteadyStateWrapper
%          extendedHovorkaModel
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

% Parameters
GEZI = p(6); % [1/min]          Glucose effectiveness
EGP0 = p(7); % [(mg/dL)/min]    Endogenous glucose production

% Meal subsystem
D1  = 0; % [g CHO]
D2  = 0; % [g CHO]

% Insulin subsystem
Isc = 0; % [mU/L] Subcutaneous insulin concentration
Ip  = 0; % [mU/L] Plasma insulin concentration

% Glucose subsystem
Ieff = 0;         % [1/min] Insulin effect
G    = EGP0/GEZI; % [mg/dL] Blood glucose concentration

% Glucose concentration measurement
Gsc = G; % [mg/dL] Subcutaneous glucose concentration

% Basal and bolus insulin flow rates
uba = 0; % [ mU/min] Insulin basal flow rate
ubo = 0; % [ mU/min] Insulin bolus flow rate

% Steady state
xs = [D1; D2; Isc; Ip; Ieff; G; Gsc];

% Steady state manipulated input
us = [uba; ubo];