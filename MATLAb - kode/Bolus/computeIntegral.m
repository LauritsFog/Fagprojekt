function i=computeIntegral(T,G)

% Compute det integral of the Glucose penalty function
% T = integration time

% This is the integal calculations
i=0;
for j=1:length(T)-1
    dt=T(j+1)-T(j);
    i = i+dt*glucosePenaltyFunction(G(j));
end

end

