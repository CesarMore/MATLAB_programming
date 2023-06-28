function dxdt = modelo3PVTOL(t,x,fq,tauT,mux)



% Ecuaciones representación de estado
m = 1;
g = 9.81;
Iyy = 0.1;
dxdt(1,1) = x(4);
dxdt(2,1) = x(5);
dxdt(3,1) = x(6);
dxdt(4,1) = (1/m)*mux;
%dxdt(4,1) = (1/m)*fq*sin(x(3));
dxdt(5,1) = -g + (1/m)*fq*cos(x(3));
dxdt(6,1) = (m/Iyy)*tauT;
