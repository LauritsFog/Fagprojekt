%% Intitialize
%  Virtual Patient (MVP) model

% Clear the command window
clc;

% Clear all variables
clear all;

% Close all figures
close all;

% Restore default settings (e.g., related to plotting)
reset(groot);

% Restore the default path (to prevent loading the wrong functions)
restoredefaultpath;

%% Load libraries
run('loadLibrary');

%% Formatting
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

%% Miscellaneous
% Conversion factors
h2min = 60;      % Convert from h   to min
min2h = 1/h2min; % Convert from min to h
U2mU  = 1e3;     % Convert from Uopen   to mU
mU2U  = 1/U2mU;  % Convert from mU  to Uopen
mmolL2mgdL = 18; % Convert from mmol/L to mg/dL

%% Create simulation scenario

% Simulation model
simModel = @mvpModel;

% Output model
outputModel = @mvpOutput;

% Observed variables
observationModel = @(t, x, p) x(7);

% Simulation method/function
simMethod = @odeEulersExplicitMethodFixedStepSize;

ctrlState = [
      0.0;  %          Initial value of integral
    108.0]; % [mg/dL] Last measurements (dummy value)

% Parameters and steady state
p = generateMVPParameters();

% Steady state time (not used)
ts = [];

% Steady state blood glucose concentration
Gs = 108; % [mg/dL]

% Compute steady state
[xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);

% If fsolve did not converge, throw an error
if(flag ~= 1), error ('fsolve did not converge!'); end

% Objective function
objectiveFunction = @asymmetricQuadraticPenaltyFunction;

% Initial and final time
t0 =  0;       % min
tf = 48*h2min; % min

% Sampling time
Ts = 5; % min

% Number of control/sampling intervals
N = (tf - t0)/Ts; % [#]

% Number of time steps in each control/sampling interval
opts.Nk = 10;

% Time span
tspan = Ts*(0:N);

% Initial condition
x0 = xs;

% Manipulated inputs
Uopen = repmat(us, 1, N);

% Disturbance variables
D = zeros(1, N);

% Disturbance variables with multiple meals
Dmealplan = MealPlan(2,false)';

% Meal and meal bolus after 1 hour
tMeal           = 1*h2min;        % [min]
idxMeal         = tMeal/Ts + 1;   % [#]
D(1, idxMeal)   = 90   /Ts;       % [g CHO/min]
Uopen(2, idxMeal)   = 0;       % [mU/min]

% Select which D to use
Duse = Dmealplan;

% Scaling factor for the objective function (purely for numerical reasons)
scalingFactor = 1e-2;

% Index of the insulin bolus in the vector of manipulated inputs
idxbo = 2;

% Initial guess of the optimal insulin bolus
ubo0 = 0; % [mU/min]

%% Simulation optimal bolus

% Controller parameters and state
ctrlPar = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.0005; %           Integral gain
      0.2500; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate 

% Halting iterations used in PID controller
haltinghours = 8;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerOptBolus;

[TOptBolus, XOptBolus, YOptBolus, UOptBolus] = closedLoopSimulationOptBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, ubo0, idxbo, scalingFactor, objectiveFunction, outputModel, opts);

% Blood glucose concentration
GscOptBolus = YOptBolus; % [mg/dL]

%% Simulation super bolus

% Controller parameters and state
ctrlPar = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.0005; %           Integral gain
      0.2500; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

[TSupBolus, XSupBolus, YSupBolus, USupBolus] = closedLoopSimulationSupBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, haltingiter, idxbo, opts);

% Blood glucose concentration
GscSupBasal = YSupBolus; % [mg/dL]

%% Visualization

c = copper(2);

% Setting critical intervals and interval colors
Gcrit = [3,3.9,10,13.9,2*13.9]*mmolL2mgdL;
Gcritcolors = {[255, 105, 105]/255;
               [255, 156, 156]/255;
               [156, 255, 159]/255;
               [255, 247, 156]/255;
               [255, 219, 156]/255};

% Create figure with absolute size for reproducibility
figure

% Plot blood glucose concentration
subplot(411);
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
plot(TOptBolus*min2h, GscOptBolus,'Color',c(1,:));
hold on
plot(TSupBolus*min2h, GscSupBasal,'Color',c(2,:));
xlim([t0, tf]*min2h);
ylim([0, max([GscOptBolus,GscSupBasal])*1.2]);
ylabel({'Blood glucose concentration', '[mg/dL]'});
% legend('Closed loop with optimal bolus')

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*Duse(1, :), 'MarkerSize', 0.1,'Color','k');
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, UOptBolus(1, [1:end, end]),'LineWidth', 4,'Color',c(1,:));
hold on
stairs(tspan*min2h, USupBolus(1, [1:end, end]),'LineWidth', 4,'Color',c(2,:));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
plot(tspan(1:end-1)*min2h, Ts*mU2U*UOptBolus(2, :),'LineStyle','none','Marker', 'o', 'MarkerSize', 4, 'Color', c(1,:));
hold on
plot(tspan(1:end-1)*min2h, Ts*mU2U*USupBolus(2, :),'LineStyle','none','Marker', '+', 'MarkerSize', 4, 'Color', c(2,:));
xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[Uopen]'});
xlabel('Time [h]');
