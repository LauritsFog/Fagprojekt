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
simModel = @mvpModel;
% simModel = @mvpNoise;

% Output model
outputModel = @mvpOutput;

% Observed variables
observationModel = @(t, x, p) x(7);

% Simulation method/function
simMethod = @odeEulersExplicitMethodFixedStepSize;
% simMethod = @odeEulerMaruyamasExplicitMethodFixedStepSize;

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
days = 10;
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

% Initial guess of the optimal insulin bolus
ubo0 = 0; % [mU/min]

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Ramping function
rampingfunction = @sigmoidRamp;

% Time before titration begins
tzero = (0.5*h2min)/Ts;

% ctrlPar = [
%       5.0;    % [min]     Sampling time
%       0.05;   %           Proportional gain
%       0.0005; %           Integral gain
%       0.2500; %           Derivative gain
%     108.0;    % [mg/dL]   Target blood glucose concentration
%     us(1)];     % [mU/min]  Nominal basal rate (overwritten below)

% Controller parameters and state
ctrlPar = [
      5.0;    % [min]     Sampling time
      0.15;   %           Proportional gain
      0.00015; %           Integral gain
      1; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate (overwritten below)

% Parameters for grid algorithm
dg = 8;
dt = 1;
gridTime = 3*h2min/Ts;

%% Simulation optimal bolus

[TOptBolus, XOptBolus, YOptBolus, UOptBolus] = closedLoopSimulationOptBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, ubo0, idxbo, scalingFactor, objectiveFunction, outputModel, rampingfunction, opts);

% Blood glucose concentration
GscOptBolus = mvpOutput(XOptBolus); % [mg/dL]

%% Simulation super bolus with PID sim.

% Computing super bolus with PID simulation
simPID = 1;

[TSupBolusPIDsim, XSupBolusPIDsim, YSupBolusPIDsim, USupBolusPIDsim] = closedLoopSimulationSupBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, tzero, ...
    haltingiter, idxbo, simPID, rampingfunction, ...
    dg, dt, gridTime, opts);

% Blood glucose concentration
GscSupBolusPIDsim = mvpOutput(XSupBolusPIDsim); % [mg/dL]

%% Simulation super bolus without PID sim.

% Computing super bolus without PID simulation
simPID = 0;

[TSupBolus, XSupBolus, YSupBolus, USupBolus] = closedLoopSimulationSupBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, tzero, ...
    haltingiter, idxbo, simPID, rampingfunction, ...
    dg, dt, gridTime, opts);

% Blood glucose concentration
GscSupBolus = mvpOutput(XSupBolus); % [mg/dL]

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
p0 = yline(ctrlPar(5),'LineWidth',1.2,'Color','r','LineStyle','--');
hold on
p1 = plot(TOptBolus*min2h, GscOptBolus,'Color',c(1,:));
hold on
p2 = plot(TSupBolus*min2h, GscSupBolus,'Color',c(2,:));
hold on
p3 = plot(TSupBolusPIDsim*min2h, GscSupBolusPIDsim,'Color',c(3,:));
xlim([t0, tf]*min2h);
ylim([0, max([GscOptBolus,GscSupBolus])*1.2]);
ylabel({'CGM measurements', '[mg/dL]'});
legend([p0,p1,p2,p3],'Steady state', 'Optimal bolus', 'Super bolus using CGM signal', 'Super bolus using PID sim.')

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*Duse(1, :), 'MarkerSize', 0.1,'Color','k');
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, UOptBolus(1, [1:end, end]),'LineWidth', 2.5,'Color',c(1,:));
hold on
stairs(tspan*min2h, USupBolus(1, [1:end, end]),'LineWidth', 2.5,'Color',c(2,:));
hold on
stairs(tspan*min2h, USupBolusPIDsim(1, [1:end, end]),'LineWidth', 2.5,'Color',c(3,:));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*UOptBolus(2, :),'filled','LineStyle','-','LineWidth', 1,'Marker', 'o', 'MarkerSize', 5, 'Color', c(1,:));
hold on
stem(tspan(1:end-1)*min2h, Ts*mU2U*USupBolus(2, :),'filled','LineStyle','-','LineWidth', 1,'Marker', 'd', 'MarkerSize', 5, 'Color', c(2,:));
hold on
stem(tspan(1:end-1)*min2h, Ts*mU2U*USupBolusPIDsim(2, :),'filled','LineStyle','-','LineWidth', 1,'Marker', '*', 'MarkerSize', 5, 'Color', c(3,:));
xlim([t0, tf]*min2h);
ylim([0, 1.2*Ts*mU2U*max([max(UOptBolus(2, :)),max(USupBolus(2, :)),max(USupBolusPIDsim(2, :))])]);
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');

%% Percent visualization

figure
subplot(1,3,1)
PlotProcent(ComputeProcent(GscOptBolus, Gcrit));
title("Optimal bolus" + newline)
subplot(1,3,2)
PlotProcent(ComputeProcent(GscSupBolus, Gcrit));
title("Super bolus with CGM slope" + newline)
subplot(1,3,3)
PlotProcent(ComputeProcent(GscSupBolusPIDsim, Gcrit));
title("Super bolus with PID sim." + newline)

%% Bolus covariance investigation

covMatrix = cov([nonzeros(UOptBolus(2,:)),nonzeros(USupBolusPIDsim(2,:)),nonzeros(USupBolus(2,:))]);
covMatrix = round(covMatrix);

figure
scatter(Ts*mU2U*nonzeros(UOptBolus(2,:)),Ts*mU2U*nonzeros(USupBolus(2,:)),'filled')
xlabel("Bolus size with optimal bolus [U]")
ylabel("Bolus size with super bolus using CGM signal [U]")

figure
scatter(Ts*mU2U*nonzeros(USupBolusPIDsim(2,:)),Ts*mU2U*nonzeros(USupBolus(2,:)),'filled')
xlabel("Bolus size with super bolus using PID sim. [U]")
ylabel("Bolus size with super bolus using CGM signal [U]")

% figure
% plot(tspan(1:100),NoiseSpikeFilter(UOptBolus(1,81:180)+randn(1,100),3))
% 
% Utemp = UOptBolus(1,81:180)+randn(1,100);
% 
% NoiseTemp = NoiseSpikeFilter(Utemp,15);
% 
% figure
% plot(tspan(1:100),lowpassfilter(NoiseTemp,5,6,100))
% 
% figure
% plot(tspan(1:100),UOptBolus(1,81:180)+randn(1,100))

%% Simulation optimal bolus standard PID

[TOptBolusStan, XOptBolusStan, YOptBolusStan, UOptBolusStan] = closedLoopSimulationOptBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, @pidController, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, ...
    ubo0, idxbo, scalingFactor, objectiveFunction, ...
    outputModel, rampingfunction, opts);

% Blood glucose concentration
GscOptBolusStan = mvpOutput(XOptBolusStan); % [mg/dL]

%% Simulation optimal bolus halting PID

[TOptBolus, XOptBolus, YOptBolus, UOptBolus] = closedLoopSimulationOptBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, @pidControllerOptBolus, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, ...
    ubo0, idxbo, scalingFactor, objectiveFunction, ...
    outputModel, rampingfunction, opts);

% Blood glucose concentration
GscOptBolus = mvpOutput(XOptBolus); % [mg/dL]

%% Simulation optimal bolus sup PID

[TOptBolusSup, XOptBolusSup, YOptBolusSup, UOptBolusSup] = closedLoopSimulationOptBolus(x0, tspan, Duse, p, ...
    simModel, observationModel, @pidControllerSupBolus, ...
    ctrlPar, ctrlState, simMethod, tzero, haltingiter, ...
    ubo0, idxbo, scalingFactor, objectiveFunction, ...
    outputModel, rampingfunction, opts);

% Blood glucose concentration
GscOptBolusSup = mvpOutput(XOptBolusSup); % [mg/dL]

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
p0 = yline(ctrlPar(5),'LineWidth',1.2,'Color','r','LineStyle','--');
hold on
p1 = plot(TOptBolusStan*min2h, GscOptBolusStan,'Color',c(1,:));
hold on
p2 = plot(TOptBolus*min2h, GscOptBolus,'Color',c(2,:));
hold on
p3 = plot(TOptBolusSup*min2h, GscOptBolusSup,'Color',c(3,:));
xlim([t0, tf]*min2h);
ylim([0, max([GscOptBolusStan,GscOptBolus])*1.2]);
ylabel({'CGM measurements', '[mg/dL]'});
legend([p0,p1,p2,p3],'Steady state', 'Optimal bolus with standard PID', 'Optimal bolus with integration halt','Optimal bolus with basal halt')

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*Duse(1, :), 'MarkerSize', 0.1,'Color','k');
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, UOptBolusStan(1, [1:end, end]),'LineWidth', 2.5,'Color',c(1,:));
hold on
stairs(tspan*min2h, UOptBolus(1, [1:end, end]),'LineWidth', 2.5,'Color',c(2,:));
hold on
stairs(tspan*min2h, UOptBolusSup(1, [1:end, end]),'LineWidth', 2.5,'Color',c(3,:));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*UOptBolusStan(2, :),'filled','LineStyle','-','LineWidth', 1,'Marker', 'o', 'MarkerSize', 5, 'Color', c(1,:));
hold on
stem(tspan(1:end-1)*min2h, Ts*mU2U*UOptBolus(2, :),'filled','LineStyle','-','LineWidth', 1,'Marker', 'd', 'MarkerSize', 5, 'Color', c(2,:));
hold on
stem(tspan(1:end-1)*min2h, Ts*mU2U*UOptBolusSup(2, :),'filled','LineStyle','-','LineWidth', 1,'Marker', '*', 'MarkerSize', 5, 'Color', c(3,:));
xlim([t0, tf]*min2h);
ylim([0, 1.2*Ts*mU2U*max([max(UOptBolusStan(2, :)),max(UOptBolus(2, :)),max(UOptBolusSup(2, :))])]);
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');

%% Percent visualization

figure
subplot(1,3,1)
PlotProcent(ComputeProcent(GscOptBolusStan, Gcrit));
title("Optimal bolus with standard PID" + newline)
subplot(1,3,2)
PlotProcent(ComputeProcent(GscOptBolus, Gcrit));
title("Optimal bolus with integration halt" + newline)
subplot(1,3,3)
PlotProcent(ComputeProcent(GscOptBolusSup, Gcrit));
title("Optimal bolus with basal halt" + newline)