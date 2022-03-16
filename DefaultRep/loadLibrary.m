function loadLibrary()
% LOADLIBRARY Add relevant library folders to the Matlab search path.
%
% SYNOPSIS:
%   AddLibrary
%
% DESCRIPTION:
%   Adds relevant library folders to the Matlab search path such that they
%   may be called from the current folder.
%
% REQUIRED PARAMETERS:
% 
% OPTIONAL PARAMETERS:
%
% RETURNS:
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

% Let the user know that the library is being loaded
fprintf('Loading diabetes library .. ');

% Add real thermodynamics functions
addpath(genpath(fullfile(pwd, './src')));

% Let the user know that the library is being loaded
fprintf('Done\n');