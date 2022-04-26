%% Bolus Calculator

%% Intro

%% sdngålsg

parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';


X_steady = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];

x0 = X_steady;
us = 25.04;

Ts = 5;
tspan = [0 Ts];

phi0 = +Inf;

U = 0:0.5:20;

D=20;

PHI = [];

for u0=U
    X=[];
    T=[];
    
    for i=1:100
        if(i==1)
            u=us+u0*1000/Ts;
            d=D;
        else
            u=us;
            d=0;
        end
        
        [ttmp,xtmp] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm);
        X = [X;xtmp];
        T = [T;ttmp];
        x0=xtmp(end,:)';
        
    end
    
    phi = computeIntegral(T,X(:,4));
    PHI = [PHI;phi];
end

figure(2)
semilogy(U,PHI)

% ?? Forstår ik hvad det her er?



