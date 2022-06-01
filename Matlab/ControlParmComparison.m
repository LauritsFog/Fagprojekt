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
simModel = @mvpModel;
% simModel = @mvpNoise;

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

% Halting iterations used in PID controller
haltinghours = 3;
haltingiter = haltinghours*h2min/Ts;

% Ramping function
rampingfunction = @sigmoidRamp;

% Time before titration begins
tzero = 0;

%% Setting parameters

% Control algorithm
ctrlAlgorithm = @pidController;

penalties = nan(5,5,5);

% Computing super bolus without PID simulation
simPID = 1;

% 0.05;   %           Proportional gain
% 0.0005; %           Integral gain
% 0.2500; %           Derivative gain

KPs = linspace(0.01, 0.3, 8);
KIs = linspace(0.00001, 0.001, 8);
KDs = linspace(0.0001, 2, 8);

%% Simulating

for i = 1:length(KPs)
    for j = 1:length(KIs)
        for k = 1:length(KDs)
            % Controller parameters and state
            ctrlPar = [
                  5.0;    % [min]     Sampling time
                  KPs(i);   %           Proportional gain
                  KIs(j); %           Integral gain
                  KDs(k); %           Derivative gain
                108.0;    % [mg/dL]   Target blood glucose concentration
                us(1)];     % [mU/min]  Nominal basal rate 
            
            [Tclosed, Xclosed, Yclosed, Uclosed] = closedLoopSimulation(x0, tspan, Duse, p, ...
                    simModel, observationModel, ctrlAlgorithm, ...
                    ctrlPar, ctrlState, simMethod, opts);

%             [TSupBolus, XSupBolus, YSupBolus, USupBolus] = closedLoopSimulationSupBolus(x0, tspan, Duse, p, ...
%                 simModel, observationModel, ctrlAlgorithm, ...
%                 ctrlPar, ctrlState, simMethod, tzero, haltingiter, idxbo, simPID, rampingfunction, opts);

            penalties(i,j,k) = asymmetricQuadraticPenaltyFunction(Tclosed,Yclosed,p);
        end
    end
end

%% Plotting

% penalties = log(penalties);
penaltiesScaled = rescale(penalties,0,255);

figure;
t = tiledlayout(2,4);
for i = 1:length(KPs)
    h(1) = nexttile;
%     Im = image([KIs(1), KDs(1)], [KIs(end), KDs(end)], reshape(penaltiesScaled(i,:,:),length(KIs),length(KDs))');
%     Im.XData = KIs;
%     Im.YData = KDs;
    image(gca, KIs, KDs, reshape(penaltiesScaled(i,:,:),length(KIs),length(KDs))');
    title("KP: " + round(KPs(i),3,"significant"))
end
cb = colorbar(h(end));
cb.TickLabels = round(linspace(min(penalties(:)), max(penalties(:)), 6), 3, "significant");
cb.Layout.Tile = 'east';
t.XLabel.String = "KI";
t.YLabel.String = "KD";
t.XLabel.FontWeight = 'bold';
t.YLabel.FontWeight = 'bold';

%% Finding optimal parameters

[v,loc] = min(penalties(:));
[ii,jj,kk] = ind2sub(size(penalties),loc);

disp(penalties(ii,jj,kk) + " at KP = " + KPs(ii) + ", KI = " + KIs(jj) + ", KD = " + KDs(kk) + " with index " + ii + ", " + kk + ", " + jj)
