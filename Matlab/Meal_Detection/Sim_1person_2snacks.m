%% Initial
clear; clc;
loadLibrary()
load('MealPlan.mat')
%Create data for one day
snack = 1;
D1 = MealPlan(1,snack);
%% Default parameters
x0=[1.24580,1.2458,0.010090,1.082110e+02,1.082110e+02,0,0];
Ts=5;
parm=[49;47;20.10;0.01060;0.00810;0.002200;1.3300;253;47;5];

%% Testing GridAlgo

%Creating the data
days=1;
% Simulation of a human in "days" days
% Need to change D1 in order to get other number of days
[t,Gm]=Simulation(x0,Ts,days,D1,parm,0);

%plot(t,Gm(:,4))
%dg=3;
%dt=1;
dg=0.1;
dt=1;
%12 is arbitary chosen the parameter doesnt do anything in the function
figure(3);
[GF,dGF,GRID]=GridAlgo(Gm(:,4),dg,dt,12,t);
plot(t,GRID*350);




