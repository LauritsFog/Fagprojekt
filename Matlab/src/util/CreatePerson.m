function p = CreatePerson()
% Description
%- This function outputs the parameters for a random person.

b = 0; % Mean

% Parameter values:

% Tau1
tau1    =  49.0;    % [min]            Insulin absorption time 
a = tau1*0.1; % Standart deviation
tau1 = tau1 + (a.*randn(1,1) + b);

tau2    =  47.0;    % [min]            Insulin absorption time
a = tau2*0.1; % Standart deviation
tau2 = tau2 + (a.*randn(1,1) + b);

p2      =   0.0106; % [1/min]          Inverse of insulin action time constant
a = p2*0.1; % Standart deviation
p2 = p2 + (a.*randn(1,1) + b);

GEZI    =   0.0022; % [1/min]          Glucose effectiveness
a = GEZI*0.1; % Standart deviation
GEZI = GEZI + (a.*randn(1,1) + b);

VG      = 253.0;    % [dL]             Glucose distribution volume
a = VG*0.1; % Standart deviation
VG = VG + (a.*randn(1,1) + b);

taum    =  47.0;    % [min]            Meal time constant
a = taum*0.1; % Standart deviation
taum = taum + (a.*randn(1,1) + b);

tausc   =   5.0;    % [min]            Subcutaneous insulin time constant
a = tausc*0.1; % Standart deviation
tausc = tausc + (a.*randn(1,1) + b);


% Parameter values held constant
CI      =  20.1;    % [dL/min]         Insulin clearance rate
SI      =   0.0081; % [(dL/mU)/min]    Insulin sensitivity
EGP0    =   1.33;   % [(mg/dL)/min]    Endogenous glucose production


% Parameter vector
p = [tau1; tau2; CI; p2; SI; GEZI; EGP0; VG; taum; tausc];
