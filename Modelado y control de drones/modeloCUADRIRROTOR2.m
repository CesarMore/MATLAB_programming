function dxdt = modeloCUADRIRROTOR2(t,x,fq,tau)

% se debe meter  en ESPACIO DE ESTADOS
m=1;
Iyy=0.1;
Ixx=0.1;
Izz=0.1;
g=9.81;
Id=[1 0 0 ; 0 1 0 ; 0 0 1];
%------------------------------------------------------------
etap = [x(10);
        x(11);
        x(12)];
 
%--------------------------------------------
% MATRIZ DE CORIOLIS
  c11=0;
 c12=(Iyy-Izz)*(x(11)*cos(x(4)) *sin(x(4)) + x(12)*sin(x(4))*sin(x(4))*cos(x(5)))+(Izz-Iyy)*x(12)*cos(x(4))* cos(x(4))*cos(x(5)) - Ixx*x(12)*cos(x(5));
 c13=(Izz-Iyy)*x(12)*cos(x(4))*sin(x(4))*cos(x(5))* cos(x(5));
 c21=(Izz-Iyy)*(x(11)*cos(x(4)) *sin(x(4)) + x(12)*sin(x(4))*sin(x(4))*cos(x(5)))+(Iyy-Izz)*x(12)*cos(x(4))* cos(x(4))*cos(x(5)) + Ixx*x(12)*cos(x(5));
 c22=(Izz-Iyy)*x(12)*cos(x(4))*sin(x(4));
 c23=-Ixx*x(12)*sin(x(5))* cos(x(5)) + Iyy*x(12)*sin(x(4))* sin(x(4))*sin(x(5))*cos(x(5)) + Izz*x(12)*cos(x(4))*cos(x(4))*sin(x(5))*cos(x(5));
 c31=(Iyy-Izz)*x(12)*sin(x(4))*cos(x(4))*cos(x(5))*cos(x(5)) - Ixx*x(11)*cos(x(5));
 c32=(Izz-Iyy)*(x(11)*sin(x(4))*cos(x(4))*sin(x(5)) + x(12)*sin(x(4))*sin(x(4))*cos(x(5)))+(Iyy-Izz)*x(12)*cos(x(4))*cos(x(4))*cos(x(5)) + Ixx*x(12)*sin(x(5))*cos(x(5))-Iyy*x(12)*sin(x(4))* sin(x(4))*sin(x(5))*cos(x(5)) - Izz*x(12)*cos(x(4))*cos(x(4))*sin(x(5))*cos(x(5));
 c33=(Iyy-Izz)*x(12)*sin(x(4))*cos(x(4))* cos(x(5))* cos(x(5))-Iyy*x(11)*sin(x(4))*sin(x(4))*sin(x(5))*cos(x(5))-Izz*x(11)*cos(x(4))*cos(x(4))*sin(x(5))*cos(x(5));
 
  C= [c11 c12 c13;
      c21 c22 c23;
      c31 c32 c33];
  
  wnt=[    1                     0                           0;
           0                cos(x(1))                -sin(x(1));
      -sin(x(2))     sin(x(1))*cos(x(2))    cos(x(1))*cos(x(2))];
 
  WN=[ 1         0            -sin(x(2));
       0    cos(x(1))   sin(x(1))*cos(x(2));
       0   -sin(x(1))   cos(x(1))*cos(x(2))];
 %------------------------ Matriz J ---------------------------------  
  j=wnt*Id*WN   ;    %3x3

%--------------------------------------------------------------------
%------------------Variables de estado -----------------------------
dxdt(1,1)=x(7);  
dxdt(2,1)=x(8);
dxdt(3,1)=x(9);
dxdt(4,1)=x(10);  
dxdt(5,1)=x(11);
dxdt(6,1)=x(12);
dxdt(7,1)=(fq/m)*cos(x(4))*sin(x(5))*cos(x(6))+sin(x(4))*sin(x(6));
dxdt(8,1)=(fq/m)*cos(x(4))*sin(x(5))*sin(x(6))-sin(x(4))*cos(x(6));
dxdt(9,1)=(fq/m)*cos(x(4))*cos(x(5))-g;
dxdt(10:12,1)=inv(j)*(tau-C*etap);
