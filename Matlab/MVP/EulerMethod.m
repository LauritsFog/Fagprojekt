function [T,X] = EulerMethod(fun,tspan,x0,N,u,d,parm)

nx = length(x0);
X = zeros(N+1,nx);
X(1,:) = x0';
x = x0;

t0 = tspan(1);
tf = tspan(2);
h = (tf-t0)/N;
T = linspace(t0,tf,N+1);

for i=1:N
    t = T(i);
    x = x+h*fun(t,x,u,d,parm);
    X(i+1,:) = x';
end

end

