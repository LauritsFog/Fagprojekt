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
run('loadLib');

%% Formatting
% Font size
fs = 15;

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
% simModel = @mvpModel;
simModel = @mvpNoise;

% Output model
outputModel = @mvpOutput;
measurementNoise = 9;

% Observed variables
observationModel = @(t, x, p) x(7);

% Simulation method/function
% simMethod = @odeEulersExplicitMethodFixedStepSize;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

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
days = 4;
hours = days*24;
t0 =  0;       % min
tf = hours*h2min; % min

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

% Disturbance variables with multiple meals
Dmealplan = MealPlan(days,true)';

% Select which D to use
Duse = Dmealplan;

% Scaling factor for the objective function (purely for numerical reasons)
scalingFactor = 1e-2;

% Index of the insulin bolus in the vector of manipulated inputs
idxbo = 2;

% Controller parameters and state
% With super bolus
% ctrlPar = [
%       5.0;    % [min]     Sampling time
%       0;   %           Proportional gain
%       0.0001; %           Integral gain
%       0.35; %           Derivative gain
%     108.0;    % [mg/dL]   Target blood glucose concentration
%     us(1)];     % [mU/min]  Nominal basal rate 

ctrlPar = [
      5.0;    % [min]     Sampling time
      0.1;   %           Proportional gain
      0.0007; %           Integral gain
      1; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate (overwritten below)

% ctrlPar = [
%       5.0;    % [min]     Sampling time
%       0.15;   %           Proportional gain
%       0.00015; %           Integral gain
%       0.5; %           Derivative gain
%     108.0;    % [mg/dL]   Target blood glucose concentration
%     us(1)];     % [mU/min]  Nominal basal rate (overwritten below)

% 0.05;   %           Proportional gain
% 0.0005; %           Integral gain
% 0.2500; %           Derivative gain

% Parameters for grid algorithm
dg = 15;
dt = 5;
gridTime = 3*h2min/Ts;
mealTime = 2.5*h2min/Ts;

% Ramping function
rampingfunction = @sigmoidRamp;

% Time before titration begins [timesteps]
tzero = (0.5*h2min)/Ts;

%% Simulation

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

[T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, mealTime, opts);

% Blood glucose concentration
Gsc = mvpOutput(X,measurementNoise); % [mg/dL]

ctrlPar = [
      5.0;    % [min]     Sampling time
      0.15;   %           Proportional gain
      0.000005; %           Integral gain
      1; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate (overwritten below)
    
% Control algorithm
ctrlAlgorithm = @pidController;

[Tclosed, Xclosed, Yclosed, Uclosed] = closedLoopSimulation(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

% Blood glucose concentration
Gscclosed = mvpOutput(Xclosed); % [mg/dL]

%% Visualization
c = copper(3);

% Setting critical intervals and interval colors
Gcrit = [3,3.9,10,13.9,2*13.9]*mmolL2mgdL;
Gcritcolors = getCritColors;

% Create figure with absolute size for reproducibility
figure

% Plot blood glucose concentration
subplot(411);
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
p2 = plot(T*min2h, Gscclosed,'Color',c(2,:));
p1 = plot(T*min2h, Gsc,'Color',c(1,:));
yline(ctrlPar(5),'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([0, max(Gsc)*1.2]);
ylabel({'CGM measurements', '[mg/dL]'});
legend([p1, p2],'ABP', 'PID pancreas')

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*Duse(1, :), 'MarkerSize', 0.1,'Color','k');
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, Uclosed(1, [1:end, end]),'LineWidth', 2.5,'Color',c(2,:));
hold on
stairs(tspan*min2h, U(1, [1:end, end]),'LineWidth', 2.5,'Color',c(1,:));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*Uclosed(2, :),'filled','LineStyle','-','LineWidth', 0.5,'Marker', 'o', 'MarkerSize', 4, 'Color', c(2,:));
hold on
stem(tspan(1:end-1)*min2h, Ts*mU2U*U(2, :),'filled','LineStyle','-','LineWidth', 0.5,'Marker', 'o', 'MarkerSize', 4, 'Color', c(1,:));
xlim([t0, tf]*min2h);
ylim([0, 1.2*Ts*mU2U*max(U(2, :))+1]);
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');

%% Percent visualization

figure

PlotProcent(ComputeProcent(Gsc, Gcrit));