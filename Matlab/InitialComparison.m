%% Initialize
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

%% Create simulation scenario
% Control algorithm
ctrlAlgorithm = @pidController;

% Simulation model
% simModel = @mvpModel;
simModel = @mvpNoise;

% Output model
outputModel = @mvpOutput;

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

% Update the nominal basal rate
ctrlPar(6) = us(1);

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

%% Simulate open loop

[Topen, Xopen] = openLoopSimulation(x0, tspan, Uopen, Duse, p, simModel, simMethod, opts);

% Blood glucose concentration
Gscopen = Xopen(7,:); % [mg/dL]

%% Simulate closed loop
% Closed-loop simulation
[Tclosed, Xclosed, Yclosed, Uclosed] = closedLoopSimulation(x0, tspan, Duse, p, ...
    simModel, observationModel, ctrlAlgorithm, ...
    ctrlPar, ctrlState, simMethod, opts);

% Blood glucose concentration
Gscclosed = Yclosed; % [mg/dL]

%% Simulate open loop with optimal bolus

for i = 1:length(tspan(1:end-1))
    if Duse(i) ~= 0 % If a meal is eaten
        % Create D array with meal at index 1 and 0's for the rest of the
        % time span
        Dtemp = zeros(1,length(tspan(i:end-1)));
        Dtemp(1) = Duse(i);
        
        % Compute the optimal bolus
        [ubo, flag] = computeOptimalBolus(ubo0, idxbo, x0, tspan(i:end-1), Uopen(:,i:end-1), Dtemp, p, ...
            scalingFactor, objectiveFunction, simModel, outputModel, simMethod, opts);

        % If fsolve did not converge, throw an error
        % if(flag ~= 1), error ('fmincon did not converge!'); end

        % Meal and meal bolus
        Uopen(idxbo, i) = ubo;
    end
end

% Simulate
[Topenbolus, Xopenbolus] = openLoopSimulation(x0, tspan, Uopen, Duse, p, simModel, simMethod, opts);

% Blood glucose concentration
Gscopenbolus = Xopenbolus(7,:); % [mg/dL]

%% Visualization

c = copper(3);

% Create figure with absolute size for reproducibility
figure;

Gcrit = [3,3.9,10,13.9,2*13.9]*18;
Gcritcolors = getCritColors;

% Plot blood glucose concentration
subplot(411);
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
yline(ctrlPar(5),'LineWidth',1.2,'Color','r','LineStyle','--');
hold on
p1 = plot(Topen*min2h,Gscopen,'Color',c(1,:));
hold on
p2 = plot(Tclosed*min2h,Gscclosed,'Color',c(2,:));
hold on
p3 = plot(Topenbolus*min2h,Gscopenbolus,'Color',c(3,:));
xlim([t0, tf]*min2h);
ylim([0, max(Gscclosed)*1.2]);
ylabel({'Blood glucose concentration', '[mg/dL]'});
legend([p1,p2,p3],'Open loop', 'Closed loop', 'Open loop with optimal bolus')

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*Duse(1, :), 'MarkerSize', 0.1,'Color','k');
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, Uopen(1, [1:end, end]),'LineWidth', 2,'Color',c(1,:));
hold on
stairs(tspan*min2h, Uclosed(1, [1:end, end]),'LineWidth', 2,'Color',c(2,:));
hold on
stairs(tspan*min2h, Uopen(1, [1:end, end]),'LineWidth', 2,'Color',c(3,:));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*Uopen(2, :),'filled','LineStyle','-','LineWidth', 0.5,'Marker', 's', 'MarkerSize', 5, 'Color', c(3,:));
xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[Uopen]'});
xlabel('Time [h]');

%% Penalty function comparison

phiopen = asymmetricQuadraticPenaltyFunction(Topen,Gscopen,p);
phiclosed = asymmetricQuadraticPenaltyFunction(Tclosed,Gscclosed,p);
phiopenoptbolus = asymmetricQuadraticPenaltyFunction(Topenbolus,Gscopenbolus,p);

figure
bar([1,2,3],[phiopen,phiclosed,phiopenoptbolus])