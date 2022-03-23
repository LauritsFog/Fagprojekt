                                    %% Simulation of MVP model

                                    
%% Intro

% Vi vil gerne simulerer hvad der sker med y(t) dvs glucose concentrationen
% i blodet når vi indskyder en given u(t) (injection af fast reacting insulin)
% i dette tilfælde er u = 25.04 mU/min
% Vi har også fået givet noget parameter parm. Disse parameter er
% tilknyttet en bestemt person og kan derfor godt være anderledes

% I denne simulering antager vi også at der ikke er givet et måltid og dermed er 
% d(t) = 0.  


%%
clear;

parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';
x0 = zeros(7,1);
fs = 14;

% Time span
tspan=[0 50000];

% 
u=25.04;
d=0;
N=50000;

%% Using an ODE solver

[T,X] = ode15s(@MVPmodel,tspan,x0,[],u,d,parm);

%% Using the Forward Euler method

[TEuler,XEuler] = EulerMethod(@MVPmodel,tspan,x0,N,u,d,parm);

%% Simulation


z = X(:,6);
Th = T/60;

ze = XEuler(:,6);
Teh = TEuler/60;

y = X(:,7);
ye = XEuler(:,7);


figure(1);

hold on 
    plot(Teh(1:1000),ze(1:1000),Th(1:75),z(1:75),"linewidth",3)
    xlabel("t [h]","fontsize",fs);
    ylabel("G [mg/dL]","fontsize",fs);
    set(gca,"fontsize",fs)
    legend("ode15s","ExplicitEulerSolver","location","northeast")
    
    plot(Teh(1:1000),ye(1:1000),Th(1:75),y(1:75),"linewidth",3)
    xlabel("t [h]","fontsize",fs);
    ylabel("y [mg/dL]","fontsize",fs);
    set(gca,"fontsize",fs)
    legend("ode15s","ExplicitEulerSolver","location","northeast")    
    
hold off

% Det kan ses at de to plottede grafer ligger oven i hinanden, hvilket
% betyder at vores EulerMethod funktion er god og brugbar.

% What is the blood glucose concentration at steady state? 
% - From the plot we see that the concentration is 110 mg/dL








