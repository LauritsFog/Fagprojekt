%% Initialization 
clear; clc; clear all;
% Function that initializes the data, and loads the functions used.
loadLib();

Days = 3;
InitData();

% 3 mmol/dl er nedre grænse
% vi vil gerne ligge mellem 3.9 og 10


%% List of functions and scripts used in the assignments

% Functions used:
%{
- loadLibrary()

- CreatePerson()

- pidController()

- mvpModel()

- mvpNoise()

- odeEulersExplicitMethodFixedStepSize()

- computeSteadyStateMVPModel()

- MealPlan(days,snack)

- closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

- ComputeProcent()

- PlotProcent()

- GridAlgo()

- GRID_Filter()

%}


%% Normal sim vs simSupBolus vs simoptbolus
% Creates the mealplan
D = MealPlan(Days,0);
D = D';

                % ---------- Normal sim --------------
% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

% Blood glucose concentration
Gsc_norm = Y; % [mg/dL]

                % ---------- optbolus sim --------------
                
% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Initial guess of the optimal insulin bolus
ubo0 = 0; % [mU/min]

% Index of the insulin bolus in the vector of manipulated inputs
idxbo = 2;

% Scaling factor for the objective function (purely for numerical reasons)
scalingFactor = 1e-2;

% Objective function
objectiveFunction = @asymmetricQuadraticPenaltyFunction;

% Output model
outputModel = @mvpOutput;

% Ramping function
rampingfunction = @sigmoidRamp;

% OptBolus calculations
%{
[TOptBolus, XOptBolus, YOptBolus, UOptBolus] = closedLoopSimulationOptBolus(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, ubo0, idxbo, scalingFactor, objectiveFunction, outputModel, rampingfunction, opts);

% Blood glucose concentration
GscOptBolus = YOptBolus; % [mg/dL]
%}                


                % ---------- SupBolus sim --------------
              
% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Computing super bolus with PID simulation
simPID = 1;

[TSupBolusPIDsim, XSupBolusPIDsim, YSupBolusPIDsim, USupBolusPIDsim] = closedLoopSimulationSupBolus(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, idxbo, simPID, rampingfunction, opts);

% Blood glucose concentration
GscSupBasalPIDsim = YSupBolusPIDsim; % [mg/dL]


                % ---------------Visualisation---------------

% Create figure with absolute size for reproducibility
figure(1);

% Plot blood glucose concentration
subplot(3,1,1)
plot(T*min2h, Gsc_norm);
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Normal simulation')

%{
subplot(3,1,2)
plot(T*min2h, GscOptBolus);
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('OptBolus simulation')
%}

subplot(3,1,3)
plot(T*min2h, GscSupBasalPIDsim);
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('SupBolus simulation')


% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

% Create the plot
figure(2);
[V_norm] = ComputeProcent(Gsc_norm, Gcrit);
subplot(1,3,1)
PlotProcent(V_norm);
title('Normal Sim')

%{
[V_opt] = ComputeProcent(GscOptBolus, Gcrit);
subplot(1,3,2)
PlotProcent(V_opt);
title('OptBolus sim')
%}

[V_sup] = ComputeProcent(GscSupBasalPIDsim, Gcrit);
subplot(1,3,3)
PlotProcent(V_sup);
title('SupBolus Sim')


%% (1) Simulate one person over 1 month with 3 meals
       
% Creates the mealplan
D = MealPlan(Days,0);
D = D';

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 0;

        % -------------------Simulation-------------------
% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,Noise); % [mg/dL]

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(1);

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

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

% Create the plot
figure(2);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);

% relevant plots not used
%{
% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*D(1, :), 'MarkerSize', 0.1);
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, U(1, [1:end, end]));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*U(2, :), 'MarkerSize', 1);
xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');
%}


%% (2) Test GRID algorithm on result from (1)

dg=20;  % 20 går rigtig godt ()
dt=20;  % 20 går rigtig godt
t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

figure(3);
plot(t,x*300)
xlim([t0, tf]*min2h);
xlabel('Time [h]');
title('Indication of when a meal is eaten')


%% Visualisation of (1) and (2)

 % -------------------Visualize-------------------
% Create figure with absolute size for reproducibility

% Create vector indicating the correct meal time spots
correctMeal = zeros(1,length(D));
for i=1:length(D)
   if D(i)>0
    correctMeal(i) = 2;
   end
end

figure(4);
subplot(2,1,1);
h = plot(t,x,'b*',t,correctMeal,'r-');
set(h,{'LineWidth'},{3;1})
xlim([t0, tf]*min2h);
ylim([0.8, 1.3])
xlabel('Time [h]');
legend('Predicted meal time','Correct meal time')
title('Compare predicted meal time and correct meal time')

subplot(2,1,2);
plot(T*min2h, Gsc,'b-',t,x*300,'r-')
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');

MealCorrectness(correctMeal,x,1)


%% (3) Simulate one person over 1 month with 3 meals and snacks
    
% Now snacks are included
snacks = 1;

% Creates the mealplan
D_snack = MealPlan(Days,snacks);
D_snack = D_snack';

        % -------------------Simulation-------------------
% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationSupBolus(x0, tspan, D_snack, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, idxbo, simPID, rampingfunction, opts);

% Blood glucose concentration
Gsc_snack = Y; % [mg/dL]

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(1);

% Plot blood glucose concentration
plot(T*min2h, Gsc_snack);
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Simulation 1 person - 31 days - 3 Meals/day')

% -------------- Evaluation of simulation -------------------

% Initialize critical range for glucose concentration in the blood    
Gcrit = [54.0000   70.2000  180.0000  250.2000  664.0593];

% Create the plot
figure(2);
[V] = ComputeProcent(Gsc_snack, Gcrit);
PlotProcent(V);


%% (4) Test GRID algorithm on result from (3)

dg=10;
dt=1;
t = T*min2h;
[GF,dGF,GRID_snack]=GridAlgo(Gsc_snack,dg,dt,12,t);
x_snack=GRID_Filter(GRID_snack);

figure(3);
plot(t,x_snack*300)
xlim([t0, tf]*min2h);
xlabel('Time [h]');
title('Indication of when a meal is eaten')


%% Visualisation of (3) and (4)

 % -------------------Visualize-------------------
% Create figure with absolute size for reproducibility

% Create vector indicating the correct meal time spots
correctMeal_snack = zeros(1,length(D_snack));
for i=1:length(D_snack)
   if D_snack(i)>0
    correctMeal_snack(i) = 2;
   end
end

figure(4);
subplot(2,1,1);
h = plot(t,x_snack,'b*',t,correctMeal_snack,'r-');
set(h,{'LineWidth'},{3;1})
xlim([t0, tf]*min2h);
ylim([0.8, 1.3])
xlabel('Time [h]');
legend('Predicted meal time','Correct meal time')
title('Compare predicted meal time and correct meal time')

subplot(2,1,2);
plot(T*min2h, Gsc_snack,'b-',t,x_snack*300,'r-')
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');

MealCorrectness(correctMeal_snack,x_snack,1)


%% (5) Simulate 100 persons over 1 month + Grid algo test

dg=10;
dt=1;
t = T*min2h;

parfor i=1:100
    
    p = CreatePerson();
    [xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);
    
    if(flag ~= 1), error ('fsolve did not converge!'); end
    % Update the nominal basal rate
    ctrlPar = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.0005; %           Integral gain
      0.2500; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    NaN];
    ctrlPar(6) = us(1);
    x0 = xs;
    
    % Creates the mealplan
    D = MealPlan(Days,0);
    D = D';
    
    % -------------------Simulation-------------------
    % Closed-loop simulation
    [T, X, Y, U] = closedLoopSimulationSupBolus(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, idxbo, simPID, rampingfunction, opts);

    % Blood glucose concentration
    Gsc = Y; % [mg/dL]
    
    [GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
    x=GRID_Filter(GRID);
    
    plot(T*min2h, Gsc,'b-',t,x*350,'r-');
    xlim([t0, tf]*min2h);
    hold on
    
end
hold off


%% (7A) Measurement noise Simulation without snack

% This is done in the function called mvpNoise

simModel = @mvpModel;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 5;

% Creates the mealplan
D = MealPlan(Days,0);
D = D';

correctMeal = zeros(1,length(D));
for i=1:length(D)
   if D(i)>0
    correctMeal(i) = 1;
   end
end

% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlParComplete, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,Noise); % [mg/dL]

% ------------------- GRID -------------------

dg= 3;
dt=10;
t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(2);
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

figure(3);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);

MealCorrectness(D,x,1)


%% (7B) Measurement noise Simulation with snack

% This is done in the function called mvpNoise

simModel = @mvpModel;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 5;

% Creates the mealplan
D = MealPlan(Days,1);
D = D';

correctMeal = zeros(1,length(D));
correctSnack = zeros(1,length(D));
for i=1:length(D)
   if D(i)>4.5
    correctMeal(i) = 300;
   end
   if D(i)<4.5 && D(i)>0
      correctSnack(i) = 300;
   end
end

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
figure(2);
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

figure(3);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);

MealCorrectness(D,x,1)


%% (8) Measurement noise simulation 100 Persons (no snack)  

simModel = @mvpModel;
% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 5;
dg=10;
dt=10;

T = 5:5:Days*24*60;
t = T*min2h;

GscAv = zeros(1,864);
GridAv = zeros(1,864);

figure(1);
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


%% (13A) EulerMaruyama noise Simulation without snack

% This is done in the function called mvpNoise

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 5;

% Creates the mealplan
D = MealPlan(Days,0);
D = D';

correctMeal = zeros(1,length(D));
for i=1:length(D)
   if D(i)>0
    correctMeal(i) = 1;
   end
end

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
figure(2);
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

figure(3);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);

MealCorrectness(D,x,1)


%% (13B) EulerMaruyama noise Simulation with snack

% This is done in the function called mvpNoise

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 5;

% Creates the mealplan
D = MealPlan(Days,1);
D = D';

correctMeal = zeros(1,length(D));
correctSnack = zeros(1,length(D));
for i=1:length(D)
   if D(i)>4.5
    correctMeal(i) = 300;
   end
   if D(i)<4.5 && D(i)>0
      correctSnack(i) = 300;
   end
end

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
figure(2);
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

figure(3);
[V] = ComputeProcent(Gsc, Gcrit);
PlotProcent(V);

MealCorrectness(D,x,1)


%% (14) EulerMaruyama noise simulation 100 Persons (no snack)  

simModel = @mvpNoise;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;
% Halting iterations used in PID controller
haltinghours = 2;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

Noise = 5;
dg=10;
dt=10;

T = 5:5:Days*24*60;
t = T*min2h;

GscAv = zeros(1,864);
GridAv = zeros(1,864);

figure(1);
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


%% Test af dg og dt parameter i Grid algo parameter:

%{
Generelt så er høje værdier af både dg og dt rigtig gode.
omkring 5-7 går det hen og bliver ligegyldigt at hæve værdierne
%}

Procent = [];
Average = [];

t1 = 15;

for i=1:t1

dg=i;  
dt=i;

t = T*min2h;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent(1,i) = x;
Average(1,i) = y;

end

figure(1);
yyaxis left
plot(1:t1,Procent)
ylabel('Procent');
yyaxis right
plot(1:t1,Average)
xlabel('dg and dt values');
ylabel('time [min]');
title('dt and dg has same value')
legend('Procent found within 1 hour','Average time taken to find meal')

% 

% --------- Test of dg---------------

% dt is now constant atthree different values:

dt1 = 5;
dt2 = 10;
dt3 = 15;
Procent1 = [];
Average1 = [];
Procent2 = [];
Average2 = [];
Procent3 = [];
Average3 = [];
t = T*min2h;
for i=1:t1

dg=i;  
dt=dt1;

[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent1(1,i) = x;
Average1(1,i) = y;

dt=dt2;

[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent2(1,i) = x;
Average2(1,i) = y;

dt=dt3;

[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent3(1,i) = x;
Average3(1,i) = y;

end


figure(2);
yyaxis left
plot(1:t1,Procent1,'r-',1:t1,Procent2,'g-',1:t1,Procent3,'b-')
ylabel('Procent');
yyaxis right
plot(1:t1,Average1,'r--',1:t1,Average2,'g--',1:t1,Average3,'b--')
xlabel('dg values');
ylabel('time [min]');
title('Test of dg')
legend('Procent found within 1 hour','Average time taken to find meal')

% Ud fra dette ses vi at blå er bedst i average min med en forskel fra ca
% 44.2 til ca. 40.5
% blå er den hvor dt = 15. så jo højere dt er jo bedre

% igen ser vi at jo højere dt og dg er jo bedre.

% --------- Test of dt---------------

% dt is now constant atthree different values:

dg1 = 5;
dg2 = 10;
dg3 = 15;
Procent1 = [];
Average1 = [];
Procent2 = [];
Average2 = [];
Procent3 = [];
Average3 = [];
t = T*min2h;
for i=1:t1

dg=dg1;  
dt=i;

[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent1(1,i) = x;
Average1(1,i) = y;

dg=dg2;

[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent2(1,i) = x;
Average2(1,i) = y;

dg=dg3;

[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);
x=GRID_Filter(GRID);
    
[x y] = MealCorrectness(D,x,1);

Procent3(1,i) = x;
Average3(1,i) = y;

end


figure(3);
yyaxis left
plot(1:t1,Procent1,'r-',1:t1,Procent2,'g-',1:t1,Procent3,'b-')
ylabel('Procent');
yyaxis right
plot(1:t1,Average1,'r--',1:t1,Average2,'g--',1:t1,Average3,'b--')
xlabel('dt values');
ylabel('time [min]');
title('Test of dt')
legend('Procent found within 1 hour','Average time taken to find meal')

% igen er blå bedst.

% alt i alt er det bedst når dg og dt er høje. 


%% test GRID Algo

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


dg=3;
dt=1;
t = T*min2h;
[GF,dGF,GRID_snack]=GridAlgo(Gsc,dg,dt,12,t);
x_snack=GRID_Filter(GRID_snack);

figure(1);
subplot(2,1,1)
plot(T*min2h, Gsc,'b-',t,x_snack*300,'r-',t,correctMeal,'g-') %,t,correctSnack,'y-')
yline(130,'LineWidth',1.2,'Color','k','LineStyle','--');
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('With Filter')

subplot(2,1,2)
plot(T*min2h, Gsc,'b-',t,GRID_snack*300,'r-',t,correctMeal,'g-') %,t,correctSnack,'y-')
yline(130,'LineWidth',1.2,'Color','k','LineStyle','--');
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Without Filter')



