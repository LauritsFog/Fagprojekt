function f = MVPmodel(t,x,in)

    u = cell2mat(in(1));
    d = cell2mat(in(2));
    parm = cell2mat(in(3));
    
    f = [d-x(1)/parm(9),
         (x(1)-x(2))/parm(9),
         u/(parm(1)*parm(3))-x(3)/parm(1),
         (x(3)-x(4))/parm(2),
         -parm(4)*x(5)+parm(4)*parm(5)*x(4),
         -(parm(6)+x(5))*x(6)+parm(7)+(1000*x(2))/(parm(8)*parm(9)),
         (x(6)-x(7))/parm(10)];
