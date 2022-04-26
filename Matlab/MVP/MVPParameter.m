function p = MVPParameter(T1,T2,Tm,Tsc,CI,p2,SI,GEZI,EGP0,VG)

% Default parameters
% T1 = 49 - Time constant [min]
% T2 = 47 - Time constant [min]
% Tm = 47 - Meal absorption time constant [min]
% Tsc = 5 - 
% CI = 20.1 - Insulin clearence [dL/min]
% p2 = 0.0106 - Inverse time constant [min^-1]
% SI = 0.0081 Insulin sensitivity [(dL/mU)/min]
% GEZI = 0.0022 - Gluco  [min^-1]
% EGP0 = 1.33  - Endogenous glucose production [mg/dL/min]
% VG = 253 - Glucose distribution volume [dL]

p = [T1, T2, CI, p2, SI, GEZI, EGP0, VG, Tm, Tsc];

end

