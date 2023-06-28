function [dxdt] = modeloRBF3(t,x,u)
% --Ejemplo1--
% f1 = x(1)^3;
% f2 = x(1)^2 + x(2)^2;
% f3 = 0;

% --Ejemplo2--
f1=2*x(1)^2*sin(x(2));
f2=x(1)^2+x(1)*x(2)+x(2)*cos(x(1));
f3=x(1)*x(3)+x(2)^2+x(3)*sin(x(2));

dxdt(1,1) = x(2) + f1;
dxdt(2,1) = x(3) + f2;
dxdt(3,1) = u + f3;

