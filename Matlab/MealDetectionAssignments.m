%% Initialization 
clear; clc;
% Function that initializes the data, and loads the functions used.
loadLibrary();

% Create a person to simulate over:
p = CreatePerson();
% -------------------Formatting----------------------
% Font size
fs = 11;

% Line width
lw = 3;

% Set default font size
set(groot, 'DefaultAxesFontSize',   fs);

% Set default line widths
set(groot, 'DefaultLineLineWidth',  lw);
set(groot, 'DefaultStairLineWidth', lw);
set(groot, 'DefaultStemLineWidth',  lw);

% Conversion factors
h2min = 60;      % Convert from h   to min
min2h = 1/h2min; % Convert from min to h
U2mU  = 1e3;     % Convert from U   to mU
mU2U  = 1/U2mU;  % Convert from mU  to U
days2h = 24;

        % -------------------Create simulation scenario-------------------
% Control algorithm
ctrlAlgorithm = @pidController;

% Simulation model
simModel = @mvpModel;

% Observed variables
observationModel = @(t, x, p) x(7);

% Simulation method/function
simMethod = @odeEulersExplicitMethodFixedStepSize;

% Controller parameters and state
ctrlPar = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.0005; %           Integral gain
      0.2500; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    NaN];     % [mU/min]  Nominal basal rate (overwritten below)
ctrlState = [
      0.0;  %          Initial value of integral
    108.0]; % [mg/dL] Last measurements (dummy value)

% Steady state time (not used)
ts = [];

% Steady state blood glucose concentration
Gs = 108; % [mg/dL]

% Compute steady state
[xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);

% If fsolve did not converge, throw an error
if(flag ~= 1), error ('fsolve did not converge!'); end

% Update the nominal basal rate
ctrlPar(6) = us(1);

% Number of days the simulation is run
Days = 3;

% Initial and final time
t0 =  0;       % min
tf = Days*days2h*h2min; % min

% Sampling time
Ts = 5; % min

% Number of control/sampling intervals
N = (tf - t0)/Ts; % [#]

% Number of time steps in each control/sampling interval
opts.Nk = 10;

% Time span
tspan = Ts*(1:N);

% Initial condition
x0 = xs;

% Statement to determin if snacks are included in the mealplan
snack = 0; % not included


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

[TOptBolus, XOptBolus, YOptBolus, UOptBolus] = closedLoopSimulationOptBolus(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, ubo0, idxbo, scalingFactor, objectiveFunction, outputModel, rampingfunction, opts);

% Blood glucose concentration
GscOptBolus = YOptBolus; % [mg/dL]
                


                % ---------- optbolus sim --------------
              
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

subplot(3,1,2)
plot(T*min2h, GscOptBolus);
xlim([t0, tf]*min2h);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('OptBolus simulation')

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

[V_opt] = ComputeProcent(GscOptBolus, Gcrit);
subplot(1,3,2)
PlotProcent(V_opt);
title('OptBolus sim')

[V_sup] = ComputeProcent(GscSupBasalPIDsim, Gcrit);
subplot(1,3,3)
PlotProcent(V_sup);
title('SupBolus Sim')


%% (1) Simulate one person over 1 month with 3 meals

        
% Creates the mealplan
D = MealPlan(Days,0);
D = D';

        % -------------------Simulation-------------------
% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

% Blood glucose concentration
Gsc = Y; % [mg/dL]

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(1);

% Plot blood glucose concentration
plot(T*min2h, Gsc);
xlim([t0, tf]*min2h);
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

dg=10;
dt=1;
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
h = plot(t,x,'b*',t,correctMeal,'r-')
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


%% (3) Simulate one person over 1 month with 3 meals and snacks
    
% Now snacks are included
snacks = 1;

% Creates the mealplan
D_snack = MealPlan(Days,snacks);
D_snack = D_snack';

        % -------------------Simulation-------------------
% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

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


%% (5) Simulate 100 persons over 1 month + Grid algo test

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
    [T, X, Y, U] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

    % Blood glucose concentration
    Gsc = Y; % [mg/dL]
    
    plot(T*min2h, Gsc);
    xlim([t0, tf]*min2h);
    hold on
    
end
hold off


%% (7) Update the MVP model with a stochastic measurement noise model

% This is done in the function called mvpNoise

% Adds the MVP model that includes the noise measurement
simModel = @mvpNoise;

% Creates the mealplan
D = MealPlan(Days,0);
D = D';

        % -------------------Simulation-------------------
% Closed-loop simulation
[T, X, Y, U] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

% Blood glucose concentration
Gsc = Y; % [mg/dL]

% -------------------Visualize-------------------
% Create figure with absolute size for reproducibility
figure(1);

% Plot blood glucose concentration
plot(T*min2h, Gsc);
xlim([t0, tf]*min2h);
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


%% (8) Simulate 100 persons over 1 month with added noise + Grid algo test

simModel = @mvpNoise;

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
    [T, X, Y, U] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

    % Blood glucose concentration
    Gsc = Y; % [mg/dL]
    
    plot(T*min2h, Gsc);
    xlim([t0, tf]*min2h);
    hold on
    
end
hold off








