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
% simModel = @mvpModel;
simModel = @mvpNoise;

% Output model
outputModel = @mvpOutput;

% Observed variables
observationModel = @(t, x, p) x(7);

% Simulation method/function
% simMethod = @odeEulersExplicitMethodFixedStepSize;
simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

ctrlState = [
      0.0;  %          Initial value of integral
    108.0]; % [mg/dL] Last measurements (dummy value)

% Steady state time (not used)
ts = [];

% Steady state blood glucose concentration
Gs = 108; % [mg/dL]

% Initial and final time
days = 31;
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

% Index of the insulin bolus in the vector of manipulated inputs
idxbo = 2;

% Parameters for grid algorithm
dg = 8;
dt = 1;
gridTime = 3*h2min/Ts;

% Ramping function
rampingfunction = @sigmoidRamp;

% Time before titration begins [timesteps]
tzero = (0.5*h2min)/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

%% Simulate

% Controller parameters and state
% Only PID
% ctrlPar = [
%       5.0;    % [min]     Sampling time
%       0.05;   %           Proportional gain
%       0.0005; %           Integral gain
%       0.25; %           Derivative gain
%     108.0;    % [mg/dL]   Target blood glucose concentration
%     us(1)];     % [mU/min]  Nominal basal rate 

% 0.05;   %           Proportional gain
% 0.0005; %           Integral gain
% 0.2500; %           Derivative gain

k = 100;

Gsc = cell(1,k);

Bolus = cell(1,k);
Basal = cell(1,k);

parfor i = 1:k
    p = CreatePerson();
    [xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);
    
    ctrlPar = [
      5.0;    % [min]     Sampling time
      0;   %           Proportional gain
      0.0001; %           Integral gain
      0.35; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate 
    
    ctrlPar(6) = us(1);
    x0 = xs;
    
    % Creates new mealplan
    D = MealPlan(days,1)';
    
    [T, X, Y, U] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, opts);

    % Blood glucose concentration
    Gsc{i} = mvpOutput(X,1)'; % [mg/dL]
    
    Bolus{i} = U(:,1)';
    Basal{i} = U(:,2)';
end

%% Percent plot

Gcrit = [3,3.9,10,13.9,2*13.9]*mmolL2mgdL;

figure
PlotProcent(ComputeProcent(cell2mat(Gsc), Gcrit));

%% Penalties

penalties = nan(1,k);
for i = 1:k
    penalties(i) = asymmetricQuadraticPenaltyFunction(tspan,Gsc{i},[]);
end

figure
histogram(penalties,20)

%% Worst case

[val, idx] = max(penalties);

figure
PlotProcent(ComputeProcent(Gsc{idx}, Gcrit));

c = copper(3);

% Setting critical intervals and interval colors
Gcritcolors = getCritColors;

% Create figure with absolute size for reproducibility
figure

% Plot blood glucose concentration
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
p1 = plot(tspan*min2h, Gsc{idx},'Color',c(1,:));
yline(Gs,'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([0, max(Gsc{idx})*1.2]);
ylabel({'CGM measurements', '[mg/dL]'});
legend([p1],'Worst case')

%% Boxplots and barplots

figure
boxplot(cell2mat(Gsc),'PlotStyle','compact')

n = 20;

figure
subplot(3,1,1)
histogram(cell2mat(Gsc),n)
title('CGM measurements')
subplot(3,1,2)
histogram(nonzeros(cell2mat(Basal)),n)
title('Basal rate')
subplot(3,1,3)
histogram(nonzeros(cell2mat(Bolus)),n)
title('Bolus size')