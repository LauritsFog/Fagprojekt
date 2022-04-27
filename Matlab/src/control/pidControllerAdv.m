function [uk, ctrlState] = pidController(tk, yk, dhatk, ctrlPar, ctrlState) %#ok
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
end