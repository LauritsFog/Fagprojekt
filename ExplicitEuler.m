function x = ExplicitEuler(f,tk,tkp1,xk,varargin)

dxk = f(tk, xk, varargin);

a = tk;
b = tkp1;

h = b-a;

x = xk + h.*dxk;

