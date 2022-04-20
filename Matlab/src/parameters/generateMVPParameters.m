function p = generateMVPParameters()
% GENERATEMVPPARAMETERS Generate parameters for the Medtronic Virtual
% Patient (MVP) model.
%
% SYNOPSIS:
%   p = generateMVPParameters()
%
% DESCRIPTION:
% A fixed set of parameters for the Medtronic Virtual Patient (MVP) model.
%
% REQUIRED PARAMETERS:
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   p - a vector of parameters (dimension: 10)
%
% p( 1); % [min]            Insulin absorption time
% p( 2); % [min]            Insulin absorption time
% p( 3); % [dL/min]         Insulin clearance rate
% p( 4); % [1/min]          Inverse of insulin action time constant
% p( 5); % [(dL/mU)/min]    Insulin sensitivity
% p( 6); % [1/min]          Glucose effectiveness
% p( 7); % [(mg/dL)/min]    Endogenous glucose production
% p( 8); % [dL]             Glucose distribution volume
% p( 9); % [min]            Meal time constant
% p(10); % [min]            Subcutaneous insulin time constant
% 
% DEPENDENCIES:
%
% See also mvpModel
%          mvpOutput
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

% Parameter values
tau1    =  49.0;    % [min]            Insulin absorption time
tau2    =  47.0;    % [min]            Insulin absorption time
CI      =  20.1;    % [dL/min]         Insulin clearance rate
p2      =   0.0106; % [1/min]          Inverse of insulin action time constant
SI      =   0.0081; % [(dL/mU)/min]    Insulin sensitivity
GEZI    =   0.0022; % [1/min]          Glucose effectiveness
EGP0    =   1.33;   % [(mg/dL)/min]    Endogenous glucose production
VG      = 253.0;    % [dL]             Glucose distribution volume
taum    =  47.0;    % [min]            Meal time constant
tausc   =   5.0;    % [min]            Subcutaneous insulin time constant

% Parameter vector
p = [tau1; tau2; CI; p2; SI; GEZI; EGP0; VG; taum; tausc];