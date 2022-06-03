%% Initialization 
clear; clc; clear all;
% Function that initializes the data, and loads the functions used.
loadLib();

Days = 3;
InitData();

% Noise Level
Noise = 5;

% Create person
p = CreatePerson();

% Creates the mealplan
Dsnack = MealPlan(Days,1);
Dsnack = Dsnack';

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
    correctMeal(i) = 300;
   end
   if Dsnack(i)<4.5 && Dsnack(i)>0
      correctSnack(i) = 300;
   end
end

% 3 mmol/dl er nedre grænse
% vi vil gerne ligge mellem 3.9 og 10

%% Measurement noise Simulation without snack

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

dg= 10;
dt=10;
t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(1);
subplot(211);
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
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 1 person - 31 days - 3 Meals/day')
hold off

subplot(212)
plot(T*min2h, Gsc,'k-',T*min2h, x*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
legend('CGM','Predicted Meal','Actual Meal')
title('Grid Algorithm on simulation')

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

figure(2);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (measurement noise - no snack)')

MealCorrectness(D,x,1)


%% Measurement noise Simulation with snack

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

dg= 10;
dt=10;
t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(4);
subplot(211);
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
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 1 person - 31 days - 3 Meals/day')
hold off

subplot(212)
plot(T*min2h, Gsc,'k-',T*min2h, x*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-',T*min2h,correctSnack*max(Gsc)*1.1,'y-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
legend('CGM','Predicted Meal','Actual Meal','Snack')
title('Grid Algorithm on simulation')

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

figure(5);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (measurement noise - snack)')

MealCorrectness(D,x,1)


%% Measurement noise simulation 100 Persons (no snack)  

simModel = @mvpModel;
% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

dg=10;
dt=10;

T = 5:5:Days*24*60;
t = T*min2h;

GscAv = zeros(1,864);
GridAv = zeros(1,864);

figure(6);
subplot(211);
hold on
for j = length(Gcrit):-1:1
        area([t0, tf]*min2h,[Gcrit(j),Gcrit(j)],'FaceColor',Gcritcolors{j},'LineStyle','none')
        hold on
end

parfor i=1:100
    
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
    x=GRID_Filter(GRID);
    
    GridAv = GridAv+x;
    
    % Plot blood glucose concentration
    plot(T*min2h, Gsc, 'Color',c(1,:)); 
    yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
    xlim([t0, tf]*min2h);
    
end
GscAv = GscAv/100;

[GF,dGF,GRID]=GridAlgo(GscAv,dg,dt,12,t);
    x=GRID_Filter(GRID);
yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');

plot(t,GscAv,'y-'); %T*min2h,GridAv*20,'r-'); %T*min2h,x*300,'r-')
xlim([t0, tf]*min2h);
ylim([0 400])
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 100 people')
hold off

subplot(2,1,2)
plot(t,GscAv,'k-',t,x*300,'r-')
xlim([t0, tf]*min2h);
ylim([min(GscAv)*0.8, max(GscAv)*1.1]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('GRID on Average of 100 people')


%% EulerMaruyama noise Simulation without snack

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

dg= 10;
dt=10;
t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(7);
subplot(211);
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
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 1 person - 31 days - 3 Meals/day')
hold off

subplot(212)
plot(T*min2h, Gsc,'k-',T*min2h, x*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
legend('CGM','Predicted Meal','Actual Meal')
title('Grid Algorithm on simulation')

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

figure(8);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (Eulermaruyama noise - no snack)')

MealCorrectness(D,x,1)


%% EulerMaruyama noise Simulation with snack

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

dg= 10;
dt=10;
t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(9);
subplot(211);
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
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 1 person - 31 days - 3 Meals/day')
hold off

subplot(212)
plot(T*min2h, Gsc,'k-',T*min2h, x*max(Gsc)*1.1,'r-',T*min2h,correctMeal*max(Gsc)*1.1,'g-',T*min2h,correctSnack*max(Gsc)*1.1,'y-')
xlim([t0, tf]*min2h);
ylim([min(Gsc)*0.8, max(Gsc)*1.1]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
legend('CGM','Predicted Meal','Actual Meal','Snack')
title('Grid Algorithm on simulation')

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

figure(10);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);
title('Preformance (Eulermaruyama noise - snack)')

MealCorrectness(D,x,1)


%% EulerMaruyama noise simulation 100 Persons (no snack)  

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;
% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

dg=10;
dt=10;

T = 5:5:Days*24*60;
t = T*min2h;

GscAv = zeros(1,864);
GridAv = zeros(1,864);

figure(11);
subplot(211);
hold on
for j = length(Gcrit):-1:1
        area([t0, tf]*min2h,[Gcrit(j),Gcrit(j)],'FaceColor',Gcritcolors{j},'LineStyle','none')
        hold on
end

parfor i=1:100
    
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
    x=GRID_Filter(GRID);
    
    GridAv = GridAv+x;
    
    % Plot blood glucose concentration
    plot(T*min2h, Gsc, 'Color',c(1,:)); 
    yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');
    xlim([t0, tf]*min2h);
    
end
GscAv = GscAv/100;

[GF,dGF,GRID]=GridAlgo(GscAv,dg,dt,12,t);
    x=GRID_Filter(GRID);
yline(ctrlParComplete(5),'LineWidth',1.2,'Color','r','LineStyle','--');

plot(t,GscAv,'y-'); %T*min2h,GridAv*20,'r-'); %T*min2h,x*300,'r-')
xlim([t0, tf]*min2h);
ylim([0 400])
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 100 people')
hold off

subplot(2,1,2)
plot(t,GscAv,'k-',t,x*300,'r-')
xlim([t0, tf]*min2h);
ylim([min(GscAv)*0.8, max(GscAv)*1.1]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('GRID on Average of 100 people')

