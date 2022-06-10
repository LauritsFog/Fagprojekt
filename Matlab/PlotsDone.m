%% Initialization 
clear; clear all;
% Function that initializes the data, and loads the functions used.
loadLib();
clc;

Days = 3;
InitData();

% Noise Level
Noise = 5;
dg= 3;
dt=1;

fs = 14;

load('TheOneAndOnlyMealPlanAndParameters.mat')

D = Dsnack;
for i=1:length(D)
   if Dsnack(i)<4.5 && Dsnack(i)>0
      D(i) = 0;
   end
end

correctMeal = zeros(1,length(D));
correctSnack = zeros(1,length(D));
for i=1:length(Dsnack)
   if Dsnack(i)>4.5
    correctMeal(i) = 1;
   end
   if Dsnack(i)<4.5 && Dsnack(i)>0
      correctSnack(i) = 1;
   end
end

% 3 mmol/dl er nedre grænse
% vi vil gerne ligge mellem 3.9 og 10


%% Measurement noise 

% This is done in the function called mvpNoise

simModel = @mvpModel;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,Noise); % [mg/dL]

% ------------------- GRID -------------------

t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
fig1 = figure(1);
fig1.Position = [50 50 1340 600];
subplot(411);
hold on
% Plot blood glucose concentration
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
plot(T*min2h, Gsc, 'Color',c(1,:)); 
yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
title('Figure 5.1','FontSize',16)
hold off

subplot(412)
plot(T*min2h, Gsc,'k-',T*min2h, GRID*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
title('Figure 5.2','FontSize',16)

%{
subplot(3,1,3)
plot(T*min2h, Gsc,'k-',t,GRID*300,'r-',t,correctMeal,'g-') %,t,correctSnack,'y-')
yline(130,'LineWidth',1.2,'Color','k','LineStyle','--');
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Without Filter','FontSize',14)
%}

% -------------- Evaluation of simulation -------------------

%{

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

fig2 = figure(2);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (measurement noise - no snack)')
%}

%fprintf('---------- Measurement noise - No snack -------------- \n \n')


% __________ Measurement noise Snack 

% This is done in the function called mvpNoise

simModel = @mvpModel;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, Dsnack, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,Noise); % [mg/dL]

% ------------------- GRID -------------------

t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility

subplot(413);
hold on
% Plot blood glucose concentration
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
plot(T*min2h, Gsc, 'Color',c(1,:)); 
yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
title('Figure 5.3','FontSize',16)
hold off

subplot(414)
plot(T*min2h, Gsc,'k-',T*min2h, GRID*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-',T*min2h,correctSnack*max(Gsc)*1.1,'y-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
legend({'CGM','Predicted Meal','Actual Meal','snack'},'Position',[0.82 0.48 0.01 0.005])
title('Figure 5.4','FontSize',16)

%{
subplot(4,1,4)
plot(T*min2h, Gsc,'b-',t,GRID*300,'r-',t,correctMeal,'g-',t,correctSnack,'y-')
yline(130,'LineWidth',1.2,'Color','k','LineStyle','--');
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Without Filter')
%}

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

%{
figure(5);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (measurement noise - snack)')
%}
%fprintf('---------- Measurement noise - With snack -------------- \n \n')

%%
saveas(fig1,[pwd '/Images/Noise.png']);


%% Measurement noise simulation 100 Persons (no snack)  

simModel = @mvpModel;
% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

T = 5:5:Days*24*60;
t = T*min2h;

GscAv = zeros(1,length(t));
GridAv = zeros(1,length(t));

M100 = 0;
N100 = 0;
Av100 = 0;

fig2 = figure(2);
fig2.Position = [50 50 1340 600];
hold on
%{
for j = length(Gcrit):-1:1
        area([t0, tf]*min2h,[Gcrit(j),Gcrit(j)],'FaceColor',Gcritcolors{j},'LineStyle','none')
        hold on
end
%}

X = 100;

parfor i=1:X
    
    %Creates new person
    p = CreatePerson();
    [xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);
    
    % Control parameters
    ctrlParComplete = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.00005; %           Integral gain
      0.5; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];
    
    ctrlParComplete(6) = us(1);
    x0 = xs;
    
    % Creates new mealplan
    D = MealPlan(Days,0);
    D = D';
    
    % -------------------Simulation-------------------
    [T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

    % Blood glucose concentration
    Gsc = mvpOutput(X,Noise); % [mg/dL]
    
    GscAv = GscAv+Gsc;
    
    [GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,T*min2h);
    x=GRID_Filter(GRID);
    
    [M, N, Av] = MealCorrectness(D,x)
    
    M100 = M100+M;
    N100 = N100+N;
    Av100 = Av100+Av;
    
    GridAv = GridAv+x;
    
    % Plot blood glucose concentration
    plot(T*min2h, Gsc, 'Color',c(1,:));
    %plot(T*min2h, GRID*400,'r-','LineWidth',0.5);
    %yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
    xlim([t0, tf]*min2h);
    
end

GscAv = GscAv/X;
M100 = M100;
N100 = N100;
procent = round(N100/M100*100,2);
Av100 = Av100/X;

A = sprintf('Total number of meals: %g',M100);
B = sprintf('Total number of founds meals: %g',N100);
B1 = sprintf('Procent meals found: %g',procent);
C = sprintf('Average time to detect meal: %g min\n',Av100);
disp(A)
disp(B)
disp(B1)
disp(C)

plot(t,GscAv,'y-'); %T*min2h,GridAv*20,'r-'); %T*min2h,x*300,'r-')
xlim([t0, tf]*min2h);
ylim([0 400])
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
%title('Simulation 100 poeple - Measurement noise')
hold off

%%
saveas(figure(2),[pwd '/Images/Noise100_31.png']);


%% EulerMaruyama noise Simulation 

% This is done in the function called mvpNoise

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,Noise); % [mg/dL]

% ------------------- GRID -------------------

t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
fig3 = figure(3);
fig3.Position = [50 50 1340 600];
subplot(411);
hold on
% Plot blood glucose concentration
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
plot(T*min2h, Gsc, 'Color',c(1,:)); 
yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
title('Figure 6.1','FontSize',16)
hold off

subplot(412)
plot(T*min2h, Gsc,'k-',T*min2h, GRID*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
title('Figure 6.2','FontSize',16)


% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

%{
figure(8);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (Eulermaruyama noise - no snack)')
%}

% :_______________EulerMaruyama noise Simulation with snack

% This is done in the function called mvpNoise

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, Dsnack, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,Noise); % [mg/dL]

% ------------------- GRID -------------------

t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
subplot(413);
hold on
% Plot blood glucose concentration
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
plot(T*min2h, Gsc, 'Color',c(1,:)); 
yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
title('Figure 6.3','FontSize',16)
hold off

subplot(414)
plot(T*min2h, Gsc,'k-',T*min2h, GRID*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-',T*min2h,correctSnack*max(Gsc)*1.1,'y-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
legend({'CGM','Predicted Meal','Actual Meal','Snack'},'Position',[0.82 0.48 0.01 0.005])
title('Figure 6.4','FontSize',16)

%%

saveas(figure(3),[pwd '/Images/EulerM.png']);

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

%{
figure(10);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (Eulermaruyama noise - snack)')

fprintf('---------- Eulermaruyama - With snack -------------- \n \n')
MealCorrectness(D,x,1)
%}


%% EulerMaruyama noise simulation 100 Persons (no snack)  

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;
% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

T = 5:5:Days*24*60;
t = T*min2h;

GscAv = zeros(1,length(t));
GridAv = zeros(1,length(t));

M100 = 0;
N100 = 0;
Av100 = 0;

fig4 = figure(4);
fig4.Position = [50 50 1340 600];
hold on
%{
for j = length(Gcrit):-1:1
        area([t0, tf]*min2h,[Gcrit(j),Gcrit(j)],'FaceColor',Gcritcolors{j},'LineStyle','none')
        hold on
end
%}

X = 100;

parfor i=1:X
    
    %Creates new person
    p = CreatePerson();
    [xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);
    
    % Control parameters
    ctrlParComplete = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.00005; %           Integral gain
      0.5; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];
    
    ctrlParComplete(6) = us(1);
    x0 = xs;
    
    % Creates new mealplan
    D = MealPlan(Days,0);
    D = D';
    
    % -------------------Simulation-------------------
    % Closed-loop simulation
    [T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

    % Blood glucose concentration
    Gsc = mvpOutput(X,Noise); % [mg/dL]
    
    GscAv = GscAv+Gsc;
    
    [GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,T*min2h);
    %x=GRID_Filter(GRID);
    
    [M, N, Av] = MealCorrectness(D,x)
    
    M100 = M100+M;
    N100 = N100+N;
    Av100 = Av100+Av;
    
    GridAv = GridAv+x;
    
    % Plot blood glucose concentration
    plot(T*min2h, Gsc, 'Color',c(1,:)); 
    %plot(T*min2h, GRID*400,'r-','LineWidth',0.5);
    %yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
    xlim([t0, tf]*min2h);
    
end
GscAv = GscAv/X;
M100 = M100;
N100 = N100;
procent = round(N100/M100*100,2);
Av100 = Av100/X;

A = sprintf('Total number of meals: %g',M100);
B = sprintf('Total number of founds meals: %g',N100);
B1 = sprintf('Procent meals found: %g',procent);
C = sprintf('Average time to detect meal: %g min\n',Av100);
disp(A)
disp(B)
disp(B1)
disp(C)

plot(t,GscAv,'y-'); %T*min2h,GridAv*20,'r-'); %T*min2h,x*300,'r-')
xlim([t0, tf]*min2h);
ylim([0 400])
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
%title('Simulation 100 people - Eulermaruyama')
hold off

%%
saveas(figure(4),[pwd '/Images/EulerM100_31.png']);

