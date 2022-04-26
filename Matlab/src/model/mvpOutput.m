function G = mvpOutput(X, p) %#ok
% MVPOUTPUT Evaluate the output, i.e., the blood glucose concentration for
% the Medtronic Virtual Patient (MVP) model.
%
% SYNOPSIS:
%   G = mvpOutput(X, p)
%
% DESCRIPTION:
% Given N+1 state vectors and the parameters in the Medtronic Virtual
% Patient (MVP) model, evaluate the blood glucose concentrations.
%
% REQUIRED PARAMETERS:
%   t - time
%   X - a vector of state variables       (dimension: 7 x N+1)
%   p - a vector parameters               (dimension: 10)
%
% OPTIONAL PARAMETERS:
%
% RETURNS:
%   G   - the blood glucose concentrations  (dimension: N+1)
%
% DEPENDENCIES:
%
% See also mvpModel
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

% Glucose subsystem
G = X(6, :); % [mg/dL] Blood glucose concentration