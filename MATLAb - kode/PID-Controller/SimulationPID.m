                        %% Simulation of PID-controller

%% Intro

% we are doing a closed loop simulation of our PID-controller




%% Initialization

% The parameters
parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';

% state vector
% angiver nogen start værdier
xs = [1.2458 1.2458 0.0101 108.2115 0 0 0]';

% Vi sætter start glucose concentrationen til 200
x0 = xs;
x0(4) = 200;

Ts = 5;
tspan = [0 Ts];

Nsteps = 500;

% Intet måltid
D=zeros(Nsteps,1);

X = [];
T = [];

%% PID controller initialization

r = 108;
y=x0(4);
y_prev = y;
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
Kp = 0.1;
Ti = 200;
Td = 10;
ii = 0;

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [ttmp,xtmp] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm);
    X = [X;xtmp];
    T = [T;ttmp];
    x0=xtmp(end,:)';
    y_prev=y;
    y=x0(6);
end

% Det kan ses fra plottene at der er intet måltid
% start koncentrationen er på 200
% Der går ca 200min / 3.3 timer indtil man er nede omkring et sikkert 
% niveau af glucose i blodet, hvilket er godt.

%% Scenario 2: Small meal

% Nu ændrer vi ikke start glucose concentrationen

x0 = xs;

% 5 min intervaller
Ts = 5;
tspan = [0 Ts];

Nsteps = 500;

% Vi tilføjer et måltid i starten på 10 g CHO /min
D=zeros(Nsteps,1);
% måltid 1
D(1) = 10;

% 8 timer senere spiser vi et måltid mere
%D(96) = 10;

X = [];
T = [];

%% PID controller initialization

r = 108;
y=x0(4);
y_prev = y;
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
Kp = 0.1;
Ti = 200;
Td = 10;
ii = 0;

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [ttmp,xtmp] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm);
    X = [X;xtmp];
    T = [T;ttmp];
    x0=xtmp(end,:)';
    y_prev=y;
    y=x0(4);
end

% Tydeligt at se at der er et måltid i starten
% Man kan se at glucose concentration starter ved en 108 og derefter
% springer op til 180 pga måltidet.
% Der går ca. 100-150 minutter til vi er nede på et normalt niveau igen,
% hvilket er hurtigt og godt...

% når man indsætter flere måltider:



%% Plot of simulation

% Plotter glucose concentrationen i blodet
figure(1);
plot(T,X(:,6))

% Plotter ???
figure(2);
plot(ttmp,xtmp(:,6))

% Plotter hvornår der bliver spist måltider
figure(3);
plot(D)
% ???

%% Scenario 3: Simulation og two days with three meals a day

% Nu ændrer vi ikke start glucose concentrationen
% Simuleringen starter kl 8
% måltider er kl 9, 13 og 18 hver dag

x0 = xs;

% 5 min intervaller
Ts = 5;
tspan = [0 Ts];

% 48 timer = 2880 min = 576 steps
Nsteps = 576;

% Måltider:
D=zeros(Nsteps,1);

% Dag 1
D(12) = 5; % kl 9
D(60) = 8; % kl 13
D(120) = 10; % kl 18

% dag 2
D(300) = 5; 
D(348) = 8;
D(408) = 10;

% måske er 10 lidt overdrevet til hvert måltid?

X = [];
T = [];

%% PID controller initialization

r = 108;
y=x0(4);
y_prev = y;
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
Kp = 0.1;
Ti = 200;
Td = 10;
ii = 0;

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [ttmp,xtmp] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm);
    X = [X;xtmp];
    T = [T;ttmp];
    x0=xtmp(end,:)';
    y_prev=y;
    y=x0(4);
end










