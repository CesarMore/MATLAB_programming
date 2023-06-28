function [dxdt] = modelo1(t,x,u)

dxdt(1,1) = x(2) + sin(x(1));
dxdt(2,1) = u; 