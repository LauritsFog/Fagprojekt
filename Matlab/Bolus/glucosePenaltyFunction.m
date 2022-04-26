function rho = glucosePenaltyFunction(G)
% Glucose penalty function

% Treshold for hypoglycemia
Gmin = 72;
% Weight parameter penalizing hypoglycemia
kappa = 1e6;
% Target glucose concentration
Gbar = 108;

% The function
rho=0.5*(G-Gbar)^2 + 0.5*kappa*max(Gmin-G,0)^2;

end

