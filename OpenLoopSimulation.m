function [Z,Y,X] = OpenLoopSimulation(f,g,h,parm,x0,T,U,D)
    N = length(T);
    M = length(x0);
    
    Z = zeros(1,N);
    Y = zeros(1,N);
    X = zeros(M,N);
    
    X(:,1) = x0;
    
    for k = 2:N-1
        X(:,k) = ExplicitEuler(f,T(k-1),T(k),X(:,k-1),U(k-1),D(k-1),parm);
        Y(1,k) = ExplicitEuler(g,T(k-1),T(k),X(7,k-1),parm);
        Z(1,k) = ExplicitEuler(h,T(k-1),T(k),X(6,k-1),parm);
    end
