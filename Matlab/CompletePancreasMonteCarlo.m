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
dg = 15;
dt = 5;
gridTime = 3*h2min/Ts;
mealTime = 2.5*h2min/Ts;

% Ramping function
rampingfunction = @sigmoidRamp;

% Time before titration begins [timesteps]
tzero = (0.5*h2min)/Ts;

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

%% Simulate

% Controller parameters and state
% 0.05;   %           Proportional gain
% 0.0005; %           Integral gain
% 0.2500; %           Derivative gain

k = 100;

GscCom = cell(1,k);
GscPID = cell(1,k);

BolusCom = cell(1,k);
BasalCom = cell(1,k);

BolusPID = cell(1,k);
BasalPID = cell(1,k);

parfor i = 1:k
    p = CreatePerson();
    [xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);

    x0 = xs;
    
    % Creates new mealplan
    D = MealPlan(days,true)';

    % For artificial pancreas complete
    ctrlPar = [
      5.0;    % [min]     Sampling time
      0.1;   %           Proportional gain
      0.0007; %           Integral gain
      1; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate 

    % Control algorithm
    ctrlAlgorithm = @pidControllerSupBolus;
    
    [TCom, XCom, YCom, UCom] = closedLoopSimulationComplete(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, idxbo, ... 
    rampingfunction, dg, dt, gridTime, mealTime, opts);
    
        % For PID controller
    ctrlPar = [
      5.0;    % [min]     Sampling time
      0.15;   %           Proportional gain
      0.000005; %           Integral gain
      1; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate (overwritten below)
    
    % Control algorithm
    ctrlAlgorithm = @pidController;

    [TPID, XPID, YPID, UPID] = closedLoopSimulation(x0, tspan, D, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

    % Blood glucose concentration
    GscCom{i} = mvpOutput(XCom,measurementNoise)'; % [mg/dL]

    BolusCom{i} = UCom(2,:)';
    BasalCom{i} = UCom(1,:)';

    % Blood glucose concentration
    GscPID{i} = mvpOutput(XPID,measurementNoise)'; % [mg/dL]
    
    BolusPID{i} = UPID(2,:)';
    BasalPID{i} = UPID(1,:)';
end

%% Penalties

penaltiesCom = nan(1,k);
for i = 1:k
    penaltiesCom(i) = asymmetricQuadraticPenaltyFunction(tspan,GscCom{i},[]);
end

penaltiesPID = nan(1,k);
for i = 1:k
    penaltiesPID(i) = asymmetricQuadraticPenaltyFunction(tspan,GscPID{i},[]);
end

n = 20;

bins = linspace(min([penaltiesCom,penaltiesPID]),max([penaltiesCom,penaltiesPID]),n);

figure
subplot(211)
histogram(penaltiesCom, bins, 'Normalization', 'probability')
xlabel("Penalty score")
ylabel("Rate [%]")
title("Penalty score distribution for complete pancreas")
subplot(212)
histogram(penaltiesPID, bins, 'Normalization', 'probability')
xlabel("Penalty score")
ylabel("Rate [%]")
title("Penalty score distribution for PID controller")

%% Percent plot

Gcrit = [3,3.9,10,13.9,2*13.9]*mmolL2mgdL;

[valCom, idxCom] = max(penaltiesCom);
[valPID, idxPID] = max(penaltiesPID);

figure
subplot(141)
PlotProcent(ComputeProcent(cell2mat(GscCom), Gcrit));
title("Complete pancreas" + newline)
subplot(142)
PlotProcent(ComputeProcent(cell2mat(GscPID), Gcrit));
title("PID controller" + newline)
subplot(143)
PlotProcent(ComputeProcent(GscCom{idxCom}, Gcrit));
title("Worst case with" + newline + "complete pancreas")
subplot(144)
PlotProcent(ComputeProcent(GscPID{idxPID}, Gcrit));
title("Worst case with" + newline + "PID controller")

%% Plotting the worst case

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
p1 = plot(tspan*min2h, GscCom{idxCom},'Color',c(1,:));
yline(Gs,'LineWidth',1.2,'Color','r','LineStyle','--');
xlim([t0, tf]*min2h);
ylim([0, max(GscCom{idxCom})*1.2]);
ylabel({'CGM measurements', '[mg/dL]'});
legend([p1],'Worst case')

%% Boxplots and barplots

% figure
% boxplot(cell2mat(Gsc),'PlotStyle','compact')

n = 20;

CGMbins = linspace(min(vec([cell2mat(GscCom);cell2mat(GscPID)])),max(vec([cell2mat(GscCom),cell2mat(GscPID)])),n);
Basalbins = linspace(min(vec([cell2mat(BasalCom);cell2mat(BasalPID)])),max(vec([cell2mat(BasalCom),cell2mat(BasalPID)])),n);

figure
subplot(3,2,1)
histogram(cell2mat(GscCom),CGMbins, 'Normalization', 'probability')
xlabel("CGM measurements [mg/dL]")
ylabel("Rate [%]")
title("Complete pancreas" + newline)
subplot(3,2,3)
histogram(cell2mat(BasalCom),Basalbins, 'Normalization', 'probability')
xlabel("Basal insulin rates [mU/min]")
ylabel("Rate [%]")
subplot(3,2,5)
histogram(Ts*mU2U*nonzeros(cell2mat(BolusCom)),n, 'Normalization', 'probability')
xlabel("Bolus insulin sizes [U]")
ylabel("Rate [%]")
subplot(3,2,2)
histogram(cell2mat(GscPID),CGMbins, 'Normalization', 'probability')
xlabel("CGM measurements [mg/dL]")
ylabel("Rate [%]")
title("PID controller" + newline)
subplot(3,2,4)
histogram(cell2mat(BasalPID),Basalbins, 'Normalization', 'probability')
xlabel("Basal insulin rates [mU/min]")
ylabel("Rate [%]")
subplot(3,2,6)
histogram(Ts*mU2U*nonzeros(cell2mat(BolusPID)),n, 'Normalization', 'probability')
xlabel("Bolus insulin sizes [U]")
ylabel("Rate [%]")

%%

totalInsulin = sum(vec(cell2mat(BasalCom))) + sum(vec(cell2mat(BolusCom)));

basalInsulinRatio = (sum(vec(cell2mat(BasalCom))))/totalInsulin;
bolusInsulinRatio = (sum(vec(cell2mat(BolusCom))))/totalInsulin;