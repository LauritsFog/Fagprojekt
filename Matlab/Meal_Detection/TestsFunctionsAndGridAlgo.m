
%% Test NoiseSpikeFilter
Gm=[1, 3, 6, 8, 11,4,30];
A=NoiseSpikeFilter(Gm,3);
disp(A)

%% Test lowPassfilter


%% Test GridAlgo
dg=3;
dt=1;

load('MealPlan.mat')
%%
x0=[1.24580,1.2458,0.010090,1.082110e+02,1.082110e+02,0,0];
Ts=5;
days=7;
parm=[49;47;20.10;0.01060;0.00810;0.002200;1.3300;253;47;5];
% Simulation af et menneske i 7 dage
[t,Gm]=Simulation(x0,Ts,days,D7,parm);
%%
%plot(t,Gm(:,4))
dg=1;
dt=1;
%12 is arbitary chosen the parameter doesnt do anything in the function
[GF,dGF,GRID]=GridAlgo(Gm(:,4),dg,dt,12,t);
figure(3)
plot(t,GRID);

