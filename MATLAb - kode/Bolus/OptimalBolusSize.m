parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';
xs = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];

x0 = xs;
us = 25.04;

Ts = 5;
tspan = [0 Ts];

phi0 = +Inf;

U = 0:0.1:15;



DD=20:2:30;

UOPT = [];

for D=DD
    D;
    phi0 = +Inf;
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

            [ttmp,xtmp] = ode45(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm);
            X = [X;xtmp];
            T = [T;ttmp];
            x0=xtmp(end,:)';

        end

        phi = computeIntegral(T,X(:,4));
        if(phi<phi0)
            phi0=phi;
            uopt=u0;
        end
    end
    UOPT = [UOPT;uopt];
end



%% q7: Least squares fitting

fs = 14;

f = figure;

% Laver et polynomium fit af 1 og anden graf til vores data
p1 = polyfit(DD'*5,UOPT,1);
p2 = polyfit(DD'*5,UOPT,2);
p3 = polyfit(DD'*5,UOPT,3);

% giver værdien af de fittede p1 og p2 til meal size stederne...
y1 = polyval(p1,DD*5);
y2 = polyval(p2,DD*5);
y3 = polyval(p3,DD*5);

% Plotter de to fittede kurver
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




