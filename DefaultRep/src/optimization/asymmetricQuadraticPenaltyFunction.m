function phi = asymmetricQuadraticPenaltyFunction(T, Z, p) %#ok
% ASYMMETRICQUADRATICPENALTYFUNCTION evaluate an objective function with an asymmetric quadratic
% penalty function.
%
% SYNOPSIS:
%   phi = asymmetricQuadraticPenaltyFunction(T, Z, p)
%
% DESCRIPTION:
% This function approximates the objective function
% 
%   phi = int_t0^tf 0.5*(G - Gbar)^2 + 0.5*kappa*max(0.0, (Gmin - G))^2 dt
% 
% using a right-rectangle rule. Here, G is the blood glucose concentration
% which depends on the outputs, Z, and the parameters, p.
%
% REQUIRED PARAMETERS:
%   T   - vector of times       (dimension: N+1)
%   Z   - matrix of outputs     (dimension: nz x N+1)
%   p   - vector of parameters  (dimension: np)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   phi - value of the objective function
%
% DEPENDENCIES:
%
% See also singleShootingObjective
%          computeOptimalBolus
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

% Setpoint for blood glucose concentration
Gbar = 108.0; % [mg/dL]

% Minimum time-in-range blood glucose concentration
Gmin = 70; % [mg/dL]

% Hypoglycemia constant
kappa = 1.0e6;

% Blood glucose concentration
G = Z(:);

% Stage cost
rho = 0.5*(G - Gbar).^2 + 0.5*kappa*max(0.0, (Gmin - G)).^2;

% Approximate the integral using a right-rectangle rule
phi = sum(rho(2:end).*diff(T(:)));