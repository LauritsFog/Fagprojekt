%% Compute the optimal bolus and simulate a single meal response for the Medtronic Virtual Patient
%  (MVP) model

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
run('../loadLibrary');

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
U2mU  = 1e3;     % Convert from U   to mU
mU2U  = 1/U2mU;  % Convert from mU  to U

%% Create simulation scenario
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

% Simulation model
simModel = @mvpModel;

% Output model
outputModel = @mvpOutput;

% Simulation method/function
simMethod = @odeEulersExplicitMethodFixedStepSize;

% Initial and final time
t0 =  0;       % min
tf = 18*h2min; % min

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
U = repmat(us, 1, N);

% Disturbance variables
D = zeros(1, N);

% Meal and meal bolus after 1 hour
tMeal           = 0*h2min;          % [min]
idxMeal         = tMeal  /Ts + 1;   % [#]
D(1, idxMeal)   = 90     /Ts;       % [g CHO/min]

% Scaling factor for the objective function (purely for numerical reasons)
scalingFactor = 1e-1;

% Index of the insulin bolus in the vector of manipulated inputs
idxbo = 2;

% Initial guess of the optimal insulin bolus
ubo0 = 0; % [mU/min]

%% Compute
% Compute the optimal bolus
[ubo, flag] = computeOptimalBolus(ubo0, idxbo, x0, tspan, U, D, p, ...
    scalingFactor, objectiveFunction, simModel, outputModel, simMethod, opts);

% If fsolve did not converge, throw an error
if(flag ~= 1), error ('fmincon did not converge!'); end

% Meal and meal bolus
U(idxbo, 1) = ubo;

% Simulate
[T, X] = openLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, opts);

% Blood glucose concentration
G = mvpOutput(X, p); % [mg/dL]

%% Visualize
% Create figure with absolute size for reproducibility
figure;

% Plot blood glucose concentration
subplot(411);
plot(T*min2h, G);
xlim([t0, tf]*min2h);
ylabel({'Blood glucose concentration', '[mg/dL]'});

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