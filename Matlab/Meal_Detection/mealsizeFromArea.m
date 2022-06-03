%% Initial
close all;clear;clc
loadLib()
load('MealPlan.mat')
%Number of days wanted in simulation
days=2;
%Create data for one day
D1 = MealPlan(days,1);
%% Default parameters
x0=[1.24580,1.2458,0.010090,1.082110e+02,1.082110e+02,0,0];
Ts=5;
parm=[49;47;20.10;0.01060;0.00810;0.002200;1.3300;253;47;5];

%% Testing GridAlgo

%Creating the data
%days=1;
% Simulation of a human in "days" days
% Need to change D1 in order to get other number of days
[t,Gm]=Simulation(x0,Ts,days,D1,parm,0);

%plot(t,Gm(:,4))
%dg=3;
%dt=1;
dg=1;
dt=1;
%12 is arbitary chosen the parameter doesnt do anything in the function
[GF,dGF,GRID]=GridAlgo(Gm(:,4),dg,dt,12,t);
figure(3)
plot(t,GRID);

%% Finding the start of the meal
meals=GRID_Filter(GRID);
figure(4)
hold on
yyaxis left
plot(t,meals,'r')
%plot(t,GRID*350,'g')
yyaxis right
plot(t,Gm(:,4),'b')
hold off


%% Meal size test
[testmeal, mealest]=MealSize(Gm(:,4),t,30);
min(testmeal)
max(testmeal)
testmeal(testmeal > 0)
plot(testmeal)
%plot(mealest)
% hold on
% yyaxis left
% plot(testmeal)
% yyaxis right
% plot(GRID)
% hold off

%% Noise
Days=2;
InitData();


[TSupBolusPIDsim, XSupBolusPIDsim, YSupBolusPIDsim, USupBolusPIDsim] = closedLoopSimulationSupBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParSupBolus, ctrlState, simMethod, tzero, haltingiter, idxbo, simPID, rampingfunction, opts);
%%
% Blood glucose concentration
GscSupBolusPIDsim = mvpOutput(XSupBolusPIDsim,4);
%GscSupBolusPIDsim = YSupBolusPIDsim; % [mg/dL]
%plot(GscSupBolusPIDsim)

[testmealNoice, mealestNoice]=MealSize(GscSupBolusPIDsim,tspan);

figure(2)
hold on
yyaxis right
plot(GscSupBolusPIDsim)
yyaxis left
plot(testmealNoice)
hold off



