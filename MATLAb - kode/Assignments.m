                            %% Part 1
                % Modeling, Simulation and Control
clear all;
clc;
addpath('PID-Controller')
addpath('MVP')
addpath('Bolus')

%% Problem 1

% x(t) = state vector 
% u(t) = injection rate of insulin [mU/min]
% d(t) = meal ingestion rate [g CHO/min]
% y(t) = measured glucose concentration [mg/dL]
% z(t) = blood glucose concentration [mg/dL]

% Parameter
parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';

% State vector
x0 = zeros(7,1);

N=50000;

% Font size
fs = 14;

% Time span
tspan=[0 50000];

% Insulin injection rate
u=25.04;

%Person is fasting
d=0;

% Solve state vector using ode15

[T,X] = ode15s(@MVPmodel,tspan,x0,[],u,d,parm);

% Blood glucose concentration is; 
z = X(:,4);

% Measured blood glucose concentration is;
y = X(:,5);

% time inteval
Th = T/60;

% Plot af både z og y
hold on 
    plot(Th(1:75),z(1:75),"linewidth",3)
    xlabel("t [h]","fontsize",fs);
    ylabel("G [mg/dL]","fontsize",fs);
    set(gca,"fontsize",fs)
    legend("ode15s","location","northeast")
    
    plot(Th(1:75),y(1:75),"linewidth",3)
    xlabel("t [h]","fontsize",fs);
    ylabel("y [mg/dL]","fontsize",fs);
    set(gca,"fontsize",fs)
    legend("ode15s","location","northeast")    
    
hold off

% The steady state vector is 

X_steady = [1.24577,1.24577,0.0101,108.211,108.211,0,0];

% steady state blood glucose:
% 108.211

% Steady state measured blood concentration:
% 108.211

%% Problem 2

% Simulate the MVP model using both ode15 and EulerMethod

%% Problem 3 - PID controller

% Implement the PIC-controller
        % PIDControl(i,r,y,y_prev,us,Kp,Ti,Td,Ts)

% Forklaring af parameter i PIDControl
%{
i = integral term, I_k
r = glucose concentration target, y-bar
y = current blood glucose concentration, y_k
y_prev = previous blood glucose concetration, y_k-1
us = insulin steady state, 
Kp Ti and Td = tuning parameters for PID
Ts = sampling time = 5 min
%}
 

                        % Closed loop simulation

% Steady state vector
X_steady = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];

% Vi sætter start glucose concentrationen til 200
x0 = X_steady;
x0(4) = 200;

% Step size
Ts = 5;
tspan = [0 Ts];

% Number of steps
Nsteps = 500;

% Intet måltid
D=zeros(Nsteps,1);

% Initialisering af T og X
X = [];
T = [];

% Insulin target
r = 108;
%current blood glucose concentration
y=x0(4);
% previous blood glucose concentration
y_prev = y;
% Insulin injection rate
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
% tuning parameters for PID
Kp = 0.1;
Ti = 200;
Td = 10;

% integral term at step i
ii = 0;

% The simulation

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
  
                        % Plot af simuleringen

% Plotter glucose concentrationen i blodet
f = figure(2);
subplot(1,2,1)
plot(T,X(:,4),ttmp,xtmp(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title('Measured blood glucose consentration',"fontsize",fs)
set(gca,"fontsize",fs)
legend("Glucose concentration","Steady state","location","northeast")

subplot(1,2,2)
plot(D,"linewidth",3)
xlabel("t [h]","fontsize",fs);
ylabel("d [g CHO/min]","fontsize",fs);
title('Meal size',"fontsize",fs)
set(f,'Position',[100 200 1100 700]);



                        % Simulation of meal size thats too big

% Start glucose concentration er nu på steady state
x0 = X_steady;

% Step size
Ts = 5;
tspan = [0 Ts];

% Number of steps
Nsteps = 500;

% Intet måltid
D=zeros(Nsteps,1);
D(1) = 10;

% Initialisering af T og X
X = [];
T = [];

% Insulin target
r = 108;
%current blood glucose concentration
y=x0(4);
% previous blood glucose concentration
y_prev = y;
% Insulin injection rate
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
% tuning parameters for PID
Kp = 0.1;
Ti = 200;
Td = 10;

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

% Plot af simuleringen

% Plotter glucose concentrationen i blodet
f = figure(3);
subplot(1,2,1)
plot(T,X(:,4),ttmp,xtmp(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title('Measured blood glucose consentration',"fontsize",fs)
set(gca,"fontsize",fs)
legend("Glucose concentration","Steady state","location","northeast")

subplot(1,2,2)
plot(D,"linewidth",3)
xlabel("t [h]","fontsize",fs);
ylabel("d [g CHO/min]","fontsize",fs);
title('Meal size',"fontsize",fs)
set(f,'Position',[100 200 1100 700]);        
        
                             %% Part 2
                    % Integration and Optimization
clear all;
clc;
addpath('PID-Controller')
addpath('MVP')
addpath('Bolus')

%% Problem 4 - Bolus Calculator

% Make a function that can compute ø = ø(x0, u0, d0, p)

% The function glucosePenaltyFunction(G) takes the glucose concentration af
% input and is used in the function computeIntegral(T,G) så compute the
% integral
    

                            % Initialisering
                            
% Parameter
parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';

% font size
fs = 14;

% Steady state vector
X_steady = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];

x0 = X_steady;

% Insulin rate
us = 25.04;

% Step size
Ts = 5;
tspan = [0 Ts];

% Bliver ikke brugt??
phi0 = +Inf;

U = 0:0.5:20;

D=20;

PHI = [];

                                   % Bolus calculations
                                   
                                   
for u0=U
    % Laver tomme vectorer
    X=[];
    T=[];
    
    for i=1:100
        if(i==1) % Laver en anden insulin rate for step 1
            u=us+u0*1000/Ts; % Ny insulin rate
            d=D; % Det bliver spist en måltid
            
            
        else % I de 99 andre step er der intet måltid og uændret insulin rate
            
            u=us;
            d=0;
            
        end
        % Løser differentialligningerne
        [ttmp,xtmp] = ode15s(@MVPModel,tspan+5*(i-1),x0,[],u,d,parm);
        X = [X;xtmp];
        T = [T;ttmp];
        x0=xtmp(end,:)';
        
    end
    % Beregner integralet af den glucose concentrations kurve som fås i
    % det inderste for-loop
    phi = computeIntegral(T,X(:,4));
    
    % Smider værdien af det integrale ind i en vector
    PHI = [PHI;phi];
end

% Plot samtlige integrale-værdier over de forskellige 
figure(2)
semilogy(U,PHI)
xlabel("Bolus size???","fontsize",fs);
ylabel("Integral-værdi","fontsize",fs);
title('Bolus calculator',"fontsize",fs)

% Det der foregår her er at man spiser et måltid på 100g CHO (dvs 20 når vi
% arbejder med 5min intevaller). Vi laver så forskellige start insulin
% rates i det første inteval og beregner integralet af den glucose
% concentration kurve fi får.
% Til sidst er integral-værdien plottet hen over de forskellige ændrede
% start insulin rates.

% af grafen vi får frem ses det at for et måltid på 100g CHO er det bedst
% at have en bolus på omkring 10-10.5


% Integrale-værdier = bolus???

%% Problem 5 - Compute the optimal bolus

% Make a function minimizing the function ø = ø(x0, u0, d0, p)
% for at given meal size

% Parameter
parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';
% Steady state vector
xs = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];
% Insulin rate
us = 25.04;
% step size
Ts = 5;
% max bolus size
umax = 15;
% meal size
D = 20;

% Function that computes the optimal bolus given the meal size
OptimalBolus(xs,us,Ts,umax,D,parm)

% This function takes the following input;
    % xs = steady state vector
    % us = insulin injection rate
    % Ts = sted size (we work with 5 normally)
    % umax = the maximal bolus size 
    % D = meal size
    % parm = parameters

% As output you get the optimal bolus size for the meal size, D, that you
% entered
    
%% Problem 6 - Optimal Bolus for different meal size

% font size
fs = 14;
% Parameter
parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';
% Steady state vector
xs = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];
% Insulin rate
us = 25.04;
% step size
Ts = 5;
% max bolus size
umax = 15;
% meal size
DD=20:2:30;

% Emty vector for the optimal boluses
UOPT = [];

% Compute optimal bolus for all the different meal sizes
for D=DD
    D;
    phi0 = +Inf;
    uopt = OptimalBolus(xs,us,Ts,umax,D,parm);
    UOPT = [UOPT;uopt];
end

% Plot the meal size and their optimal bolus
plot(DD,UOPT,"linewidth",3)
xlabel("Meal Size","fontsize",fs);
ylabel("Optimal Bolus","fontsize",fs);
title('Optimal Bolus for differnt meal size',"fontsize",fs)

%% Problem 7 - Least square fit    


f = figure;
% Laver et polynomium fit af 1,2 og 3 grad af vores data
p1 = polyfit(DD'*5,UOPT,1);
p2 = polyfit(DD'*5,UOPT,2);
p3 = polyfit(DD'*5,UOPT,3);

% giver værdien af de fittede p1, p2 og p3 til meal size stederne...
y1 = polyval(p1,DD*5);
y2 = polyval(p2,DD*5);
y3 = polyval(p3,DD*5);

% Plotter de tre fittede grafer sammen med den "korrekte" graf
subplot(1,3,1)
plot(DD'*5,UOPT,DD'*5,y1)
title('linear fit (a*x + b)',"fontsize",fs)
xlabel("Meal size","fontsize",fs);
ylabel("Optimal Bolus","fontsize",fs);

subplot(1,3,2)
plot(DD'*5,UOPT,DD'*5,y2)
title('Quadratic fit (a*x^2 + b*x + c)',"fontsize",fs)
xlabel("Meal size","fontsize",fs);
ylabel("Optimal Bolus","fontsize",fs);

subplot(1,3,3)
plot(DD'*5,UOPT,DD'*5,y3)
title('Cubic fit (a*x^3 + b*x^2 + c*x + d)',"fontsize",fs)
xlabel("Meal size","fontsize",fs);
ylabel("Optimal Bolus","fontsize",fs);
set(f,'Position',[100 200 1100 700]);


% Af grafen ses det at det lineare fit ikke er helt godt. anden og tredje
% grad fit er meget fint. selvfølgelig er tredje grad bedere, men anden
% graf er skam også udemærket.


                             %% Part 3
                         % Test of Parameters

clear all;
clc;
                         
%% Initialisering                        

% Default parameter
parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';

% Start glucose concentration er nu på steady state
x0 = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];
xs = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];
% Step size
Ts = 5;
tspan = [0 Ts];

% font size
fs = 14;

% Number of steps
Nsteps = 500;

% Intet måltid
D=zeros(Nsteps,1);
D(1) = 10;

% Initialisering af T og X
X = [];
T = [];

% Insulin target
r = 108;
%current blood glucose concentration
y=x0(4);
% previous blood glucose concentration
y_prev = y;
% Insulin injection rate
us = 25.04;
%us = 20; %You can also try the case where the basal insulin infusion rate is incorrect
% tuning parameters for PID
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

% Plot af simuleringen

% Plotter glucose concentrationen i blodet

figure(1)
plot(T,X(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title('Measured blood glucose consentration',"fontsize",fs)
set(gca,"fontsize",fs)


% for bolus calculations
umax = 15;

%% Test af forskellige parameter

% vi starter med at lave en function der indeholder parameterne, så man
% nemt kan skifte dem ud.

% MVPParameters som ligger inde under MVP mappen gør dette

% Input er de 10 forskellige parameter
% MVPParameter(T1.T2,Tm,Tsc,CI,p2,SI,GEZI,EGP0,VG)

% Output er en vektor med disse parameter

T1 = 49;
T2 = 47;
Tm = 47; 
Tsc = 5; 
CI = 20.1;
p2 = 0.0106; 
SI = 0.0081; 
GEZI = 0.0022; 
EGP0 = 1.33;
VG = 253;

% Hvad er ændret???????
string1 = 'VG Fordoblet';
string2 = 'VG Halveret';

parm1 = MVPParameter(T1,T2,Tm,Tsc,CI,p2,SI,GEZI,EGP0,VG*2);
parm2 = MVPParameter(T1,T2,Tm,Tsc,CI,p2,SI,GEZI,EGP0,VG/2);

% Initialisering af T og X
Xnew = [];
Tnew = [];
Xnew1 = [];
Tnew1 = [];
ii = 0;

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [tnew,xnew] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm1);
    Xnew = [Xnew;xnew];
    Tnew = [Tnew;tnew];
    x0=xnew(end,:)';
    y_prev=y;
    y=x0(4);
end

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [tnew1,xnew1] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm2);
    Xnew1 = [Xnew1;xnew1];
    Tnew1 = [Tnew1;tnew1];
    x0=xnew1(end,:)';
    y_prev=y;
    y=x0(4);
end

f = figure(2);
subplot(3,1,1)
plot(T,X(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title('Default Parameter',"fontsize",fs)

subplot(3,1,2)
plot(Tnew,Xnew(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title(string1,"fontsize",fs)

subplot(3,1,3)
plot(Tnew1,Xnew1(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title(string2,"fontsize",fs)
set(f,'Position',[1700 100 600 1500]); % [x, y , bredte, højde]

% Undersøgelse af bolus
% Meal size = 20
% Max bolus
umax = 20;
disp(string1)
OptimalBolus(xs,us,Ts,umax,20,parm1)

disp(string2)
OptimalBolus(xs,us,Ts,umax,20,parm2)

%% Noter til undersøgelse af parameter

% Default parameter info
%{
 Optimal bolus ved default parameter: 9.40
Glucose concentration til forskellige tidspunkter for default parameter
t(1000) = 101.181
t(1500) = 105.226
t(2000) = 106.929

T1 = 49 - Time constant [min]
T2 = 47 - Time constant [min]
Tm = 47 - Meal absorption time constant [min]
Tsc = 5 - ?????
CI = 20.1 - Insulin clearence [dL/min]
p2 = 0.0106 - Inverse time constant [min^-1]
SI = 0.0081 Insulin sensitivity [(dL/mU)/min]
GEZI = 0.0022 - Gluco  [min^-1]
EGP0 = 1.33  - Endogenous glucose production [mg/dL/min]
VG = 253 - Glucose distribution volume [dL]
%}

% Undersøgelse af de 10 forskellige parameter
%{
                    T1 Undersøgelse:
Fordobling:
- Ingen ændring i max peak af G
- Det tager længere tid at ramme 108 mg/dL
- den dykker lidt længere ned (rammer næsten 90)
- opt bolus = 9.6

Halvering
- Ingen ændring i max peak af G
- det tager kortere tid at ramme 108
- dykker ikke helt ned til 90 som fordoblingen
- optimal bolus = 9.2

                    T2 Undersøgelse:
Fordobling:
- Tager også lidt længere tid, men ikke meget
- optimal bolus 9.6

Halvering
- Tager lidt kortere tid, men ikke meget
- optimal bolus = 9.2

- Ingen ændring i max peak af G

                    Tm Undersøgelse:

Fordobling:
- Max peak 154-155
- optimal bolus 8.9
- Gør også at det tager længere tid at ramme 108

Halvering
- Max peak 210 ish
- optimal bolus = 6.7
- det tager kortere tid at ramme 108

- Basically ændrer den her mest på max peak

                    Tsc Undersøgelse:

Fordobling:
- Ingen ændring i noget

Halvering
- Ingen ændring i noget

- Denne parameter bliver ikke brugt til noget???

                    
                        CI Undersøgelse:

Fordobling:
- max peak 190 
- Den kommer aldrig rigtig ned på 108
- Det tager i hvert fald voldsomt lang tid
- Optimal bolus = returnerer max bolus, så den er meget høj

Halvering
- max peak 160
- Dykker helt ned til omkring 40, hvilket ikke er godt
- tar omkring 2000 min at komme op på 108, hvilket ikke er godt
- Den har også et ret kraftigt fald hvilket nok ikke er særlig rart for
mennekser at opleve.
- Optimal bolus = returnerer 0 ???

- Generelt tror jeg ikke denne parameter skal påvirket.
- Den er også "Insulin clearence" så tror den er ret fastsat


                    p2 Undersøgelse:
Fordobling:
- generelt en smule hurtigere i starten, men langsommere ved t = 2000
t(1000) = 101.856
t(1500) = 105.313
t(2000) = 106.851
- Optimal bolus = 8.8

Halvering
- generelt langsommere i starten, men hurtigere ved t= 2000
t(1000) = 98.52
t(1500) = 104.97
t(2000) = 107.238
- Optimal bolus = 10

Generelt
- Intet ændring i max peak
- ligner meget default parameter grafen
- antager at der ikke er den store signifikante ændring


                           SI Undersøgelse:

Fordobling:
- er stort set lige så hurtig som default til at ramme 108
t(1000) = 101.326
t(1500) = 106.725
t(2000) = 107.772
- Max peak på 123-124 ish
- dykker helt ned til omkring 60, (ikke godt)
        - indikerer at U i starten er voldsomt høj
- Optimal bolus = 0

Halvering
- Meget langsommere om at ramme 108. rammer faktisk ikke
t(1000) = 125.418
t(1500) = 117.615
t(2000) = 113.61
- max peak næsten på 250. får aldrig et dyk
            - kan indikerer at u aldrig bliver særlig stor
- Optimal bolus = 20

Generelt
- stor forskel
- virker lidt som om den påvirker hvor stort u bliver under simulationen
- SI = Insulin sensitivity
- giver røv god mening at når man så fordobler den virker den voldsomt godt
og når man halverer den virker den nærmest ikke

- Tror også dette er en af dem der er fastlagt og ikke skal ændres

                        GEZI Undersøgelse:

Fordobling:
en anelse langsommere end default, men ikke af betydning (tænker jeg)
t(1000) = 98.657
t(1500) = 103.911
t(2000) = 106.39
- max peak 155 ish
- dykker ned til omkring 86-87
- Optimal bolus = 7.4

Halvering
- hurtigere end default
t(1000) = 105.06
t(1500) = 106.649
t(2000) = 107.549
- max peak næsten på 220
- rammer 100 ved dyk
- Optimal bolus = 10.6

Generelt
- ved ikke helt hvad denne parameter beskriver
- Den påvirker i hvert fald max peak og hvor meget den dykker ned efter max
peak
- Generelt ligner grafen ret meget default grafen


                        EGP0 Undersøgelse:
- Endogenous glucose production

Fordobling:
meget langsommere. rammer ikke helt 108 til en god tid
t(1000) = 115.671
t(1500) = 112.633
t(2000) = 110.877
- max peak er på næsten 350...
- optimal bolus = 20

Halvering
- den her er helt fucked. Den dykker ned på 35 med det samme
- har intet peak
t(1000) = 52.31
t(1500) = 79.98
t(2000) = 103.627
- Optimal bolus = 0

Generelt
- Den her er helt fucked...
- tror ikke den her parameter skal fuckes med...
- den er også fastsat forhåbentlig

                        VG Undersøgelse:
-  Glucose distribution volume

Fordobling:
- En anelse hurtigere om at ramme 108
t(1000) = 104.891
t(1500) = 106.755
t(2000) = 107.519
- max peak er på næsten 145
- min er på omkring 100
- optimal bolus = 4.6

Halvering
- langsommere om at ramme 108
t(1000) = 94.38
t(1500) = 101.815
t(2000) = 105.44
- max peak er på omkring 250
- min er på omkring 80
- Optimal bolus = 14.5

Generelt
- Den påvirker max peak
- det virker som om den også påvirker hvor meget der dykkes ned
%}

% Hvad har vi lært:
%{

- De to tidskonstanter T1 og T2 ændrer ikke så meget på bolus og ændrer en
smule på hastigheden af stabilisering af G (hvor hurtigt den rammer 108)
De ændrer overhoved ikke på max peak

- Tm, Meal absorption, ændrer på max peak mest. Selve grafen for G ligner
meget den for default. Hastighed for stabilitet er også ændret en smule.
bolus går op når Tm går op og omvendt

- CI, insulin clearence, har voldsomme udsving når man ændrer den. Det kan
godt tyde på at dette er en parameter der er fastlagt og ikke skal ændres
Den skal i hvert fald ikke fordobles eller halveres

- p2, Inverse time constant, Ændrer ikke på max peak og ændrer en smule på
hastigheden af stabilitet. dog ikke signifikant
Den er lidt mærkelig, da der ikke sker det store ved at halverer eller
fordoble den.
fordobling = lavere bolus
halvering = højere bolus

- SI, insulin sensitivity, ændrer rigtig meget på max peak og hvor meget
grafen dykker efter peak.
virker lidt som om den påvirker hvor stort u bliver under simulationen
giver røv god mening at når man så fordobler den virker den voldsomt godt
og når man halverer den virker den nærmest ikke
Tror også dette er en af dem der er fastlagt og ikke skal ændres
bolus er helt fucked på den her

- GEZI påvirker også max peakog hvor meget den dykker ned efter
bolus bliver mindre når den fordobles og større når den halveres
Generelt ligner grafen ret meget default grafen

- EGP0. Den her er helt fucked...
Halvering = dykker ned på 35 ned det samme
tror ikke den her parameter skal fuckes med...
den er også fastsat forhåbentlig

VG, Glucose distribution volume, ændrer bolus en del
fordobling = meget lav bolus, halvering = høj bolus
Den påvirker max peak
det virker som om den også påvirker hvor meget der dykkes ned

%}

%% Forsøg - forbedre stabilisering tiden

p = 0.5;
% Parameter værdier
T1 = 49 -49*p;
T2 = 47 -47*p;
Tm = 47 -47*p; 
Tsc = 5; 
CI = 20.1;
p2 = 0.0106 + 0.0106*p; 
SI = 0.0081; 
GEZI = 0.0022 - 0.0022*p; 
EGP0 = 1.33;
VG = 253 + 253*p;

% Hvad er ændret???????
string3 = 'T1,T2,Tm = ned VG = op';

% parameter
parm3 = MVPParameter(T1,T2,Tm,Tsc,CI,p2,SI,GEZI,EGP0,VG);

% Initialisering af T og X
Xopt = [];
Topt = [];

ii = 0;

for i=1:Nsteps
    [u,ii] = PIDControl(ii,r,y,y_prev,us,Kp,Ti,Td,Ts);
    u = max(u,0);
    d = D(i);
    [topt,xopt] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm3);
    Xopt = [Xopt;xopt];
    Topt = [Topt;topt];
    x0=xopt(end,:)';
    y_prev=y;
    y=x0(4);
end

f = figure(2);
subplot(2,1,1)
plot(T,X(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title('Default Parameter',"fontsize",fs)

subplot(2,1,2)
plot(Topt,Xopt(:,4),"linewidth",3)
xlabel("t [min]","fontsize",fs);
ylabel("G [mg/dL]","fontsize",fs);
title(string3,"fontsize",fs)
set(f,'Position',[1700 100 600 1500]); % [x, y , bredte, højde]

%{
10% ændring
T1,T2,Tm = ned VG = op'

t(1000) = 102.069
t(1500) = 105.585
t(2000) = 107.04

20% ændring
T1,T2,Tm = ned VG = op'

t(1000) = 102.816
t(1500) = 105.84
t(2000) = 107.13

30% ændring
T1,T2,Tm = ned VG = op'

t(1000) = 103.351
t(1500) = 106.045
t(2000) = 107.212

30% ændring
T1,T2,Tm = ned VG = op'

t(1000) = 103.83
t(1500) = 106.215
t(2000) = 107.261

- Go længere man får disse værdier ned / op jo bedre bliver tiden.
- peak falder en smule 
______________________________________________________________________

ændring 40%
T1,T2,Tm, GEZI = ned VG, p2 = op'

t(1000) = 105.443
t(1500) = 106.871
t(2000) = 107.525

ændring 50%
T1,T2,Tm, GEZI = ned VG, p2 = op'

t(1000) = 106.625
t(1500) = 107.392
t(2000) = 107.731

- peak lidt lavere end 180
%}

















