%% Initial
loadLibrary()
load('MealPlan.mat')
%Create data for one day
D1 = MealPlan(1);
%% Default parameters
x0=[1.24580,1.2458,0.010090,1.082110e+02,1.082110e+02,0,0];
Ts=5;
parm=[49;47;20.10;0.01060;0.00810;0.002200;1.3300;253;47;5];




%% Testing GridAlgo

%Creating the data
days=1;
% Simulation of a human in "days" days
% Need to change D1 in order to get other number of days
[t,Gm]=Simulation(x0,Ts,days,D1,parm,1);

figure(1)
plot(t,Gm(:,4),"linewidth",3)
xlabel("t [min]");
ylabel("G [mg/dL]");
title('Measured blood glucose consentration')
set(gca)
legend("Glucose concentration","Steady state","location","northeast")

%plot(t,Gm(:,4))
%dg=3;
%dt=1;
dg=1;
dt=1;
%12 is arbitary chosen the parameter doesnt do anything in the function
[GF,dGF,GRID]=GridAlgo(Gm(:,4),dg,dt,12,t);
plot(t,GRID*350);


