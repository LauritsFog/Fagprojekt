parm = [49,47,20.1,0.0106,0.0081,0.0022,1.33,253,47,5];

tspan = linspace(1,24*60,24*60);
D = zeros(1,length(tspan));
U = ones(1,length(tspan)).*25.04;

x0 = zeros(1,7);

[Z,Y,X] = OpenLoopSimulation(@MVPmodel,@CGMsensor,@GluCon,parm,x0,tspan,U,D);

% X(6,:) = [mg/dL] is the blood glucose concentration in the blood. 
% X(7,:) = [mg/dL] is the measured subcutaneous glucose concentration.

%%

plot(tspan,X(6,:),'-');

%%

[T,X] = ode15s(@MVPmodel,[0 24*60],x0',[U(1),D(1),parm]);

%%

