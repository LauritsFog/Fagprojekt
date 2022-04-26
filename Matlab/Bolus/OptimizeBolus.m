parm = [49 47 20.1 0.0106 0.0081 0.0022 1.33 253 47 5]';

xs = [1.2458, 1.2458 , 0.01009, 108.211, 108.211, 0, 0 ];

x0 = xs;
us = 25.04;

Ts = 5;
tspan = [0 Ts];

phi0 = +Inf;

U = 0:0.1:15;



D=20;

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
    if(phi<phi0)
        phi0=phi;
        uopt=u0;
    end
end

disp('Optimal bolus (in U):')
uopt