function uopt = OptimalBolus(xs,us,Ts,umax,D,parm)

% Function for finding the bolus size that minimises ø
% xs = steady state vector
% us = insulin injection rate
% Ts = sted size (we work with 5 normally)
% umax = the maximal bolus size 
% D = meal size
% parm = parameters


x0 = xs;

tspan = [0 Ts];

phi0 = +Inf;

U = 0:0.1:umax;

for u0=U
    % creates emty vectors for X and T
    X=[];
    T=[];
    
    for i=1:100
        if(i==1) % updates the insulin rate when we eat our meal
            u=us+u0*1000/Ts;
            d=D;
        else
            u=us;
            d=0;
        end
        
        % Solves the differential equations and computes glucose
        % concentration curve
        [ttmp,xtmp] = ode15s(@MVPmodel,tspan+5*(i-1),x0,[],u,d,parm);
        X = [X;xtmp];
        T = [T;ttmp];
        x0=xtmp(end,:)';
        
    end
    
    % Computes the integral of the glucose concentration curve
    phi = computeIntegral(T,X(:,4));
    
    % Updates the optimal bolus size as the foor-loop progresses
    if(phi<phi0)
        phi0=phi;
        uopt=u0;
    end
end

end

