function [T , X] = Simulation(x0,Ts,days,D,parm,plots)

%{
x0 = steady state vector
Ts = step size (minutes)
nd = number of days we simulate
r = insulin target
us = insulin injection rate
D = Meal vector
parm = parameters
%}

% Time steps
tspan = [0 Ts];

Nsteps = days*24*60/Ts;

% Initialisering af T og X
X = [];
T = [];

% Insulin target
r = 108;
%current blood glucose concentration
y=x0(4);
% previous blood glucose concentration
y_prev = y;
% Insulin injection rate
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
% tuning parameters for PID
Kp = 0.1;
Ti = 200;
Td = 10;

% integral term at step i
ii = 0;

% The simulation

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [ttmp,xtmp] = ode15s(@MVPModel,tspan+5*(i-1),x0,[],u,d,parm);
    X = [X;xtmp];
    T = [T;ttmp];
    x0=xtmp(end,:)';
    y_prev=y;
    y=x0(4);
end

%% Parameter deciding to plot or not. Leaving plots blank gives the plot
if ~exist('plots','var')
     % third parameter does not exist, so default it to something
      plots = 0;
 end



if plots==1

else
    % Plotter glucose concentrationen i blodet
    fs = 14;
                

figure(1)
plot(T,X(:,4), ttmp,xtmp(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title('Measured blood glucose consentration',"fontsize",fs)
set(gca,"fontsize",fs)
legend("Glucose concentration","Steady state","location","northeast")

figure(2)
plot(D,"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("d [g CHO/min]","fontsize",fs);
title('Meal size',"fontsize",fs)
%%set(f,'Position',[100 200 1100 700]);
end



end

