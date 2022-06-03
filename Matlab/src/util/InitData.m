
%INITDATA Summary of this function goes here
%   Detailed explanation goes here


        % --------------------- Formatting --------------------- 

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


        % --------------------- Miscellaneous --------------------- 

% Conversion factors
h2min = 60;      % Convert from h   to min
min2h = 1/h2min; % Convert from min to h
U2mU  = 1e3;     % Convert from Uopen   to mU
mU2U  = 1/U2mU;  % Convert from mU  to Uopen
mmolL2mgdL = 18; % Convert from mmol/L to mg/dL
days2h = 24;

 % --------------------- Create simulation scenario --------------------- 


% Simulation model
%simModel = @mvpModel;
 simModel = @mvpNoise;

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

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

% Control algorithm
ctrlAlgorithm = @pidControllerSupBolus;

% Ramping function
rampingfunction = @sigmoidRamp;

% Time before titration begins
tzero = 0;

% Computing super bolus with PID simulation
simPID = 1;

 % --------------------- Control Parameters --------------------- 

% Controller parameters and state
ctrlParOptBolus = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.0005; %           Integral gain
      0.2500; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate 

% Controller parameters and state
ctrlParSupBolus = [
      5.0;    % [min]     Sampling time
      0.43;   %           Proportional gain
      1e-05;  %           Integral gain
      2;      %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    us(1)];     % [mU/min]  Nominal basal rate 



