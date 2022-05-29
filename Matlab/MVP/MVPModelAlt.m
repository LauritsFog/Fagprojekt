function dx = MVPModelAlt(t,x,u,d,parm)

tau1 = parm(1);
tau2 = parm(2);
CI = parm(3);
p2 = parm(4);
SI = parm(5);
GEZI = parm(6);
EGP0 = parm(7);
Vg = parm(8);
taum = parm(9);
tsc = parm(10); %%%%%


dx = zeros(6,1);
dx(1) = u/(tau1*CI)-x(1)/tau1;
dx(2) = (x(1)-x(2))/tau2;
dx(3) = -p2*x(3)+SI*p2*x(2);
dx(4) = -(GEZI+x(3))*x(4)+EGP0+1000*x(7)/(Vg*taum);
dx(5) = (x(4)-x(5))/tsc; %%%%%
dx(6) = d-x(6)/taum;
dx(7) = (x(6)-x(7))/taum;

