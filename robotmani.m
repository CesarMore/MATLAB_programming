function [Xp,u]=robotmani(t,X)

%Variables de estado
x1 = X(1);
x2 = X(2);
x3 = X(3);
x4 = X(4);
%vector de estados
x = [x1;x2;x3;x4];
%Variables de estado estimados
xh1 = X(5);
xh2 = X(6);
xh3 = X(7);
xh4 = X(8);
%vector de estados estiamados
xh = [xh1;xh2;xh3;xh4];

%Parametros de simulacipon 
mgL = 10;
I = 1;
J = 1;
k = 100;
a = mgL/I;
b = k/I;
c = k/J;
d = 1/J;

%Cambio de variables
z1 = x1;
z2 = x2;
z3 = -a*sin(x1)-b*(x1-x3);
z4 = -a*x2*cos(x1)-b*(x2-x4);
%vector de cambio de variables 
z = [z1;z2;z3;z4];

%Ganacias 
k1 = 0.0937;
k2 = 1.7969;
k3 = 5;
k4 = 4.0625;
% k1 = 4788;
% k2 = 2319;
% k3 = 420;
% k4 = 34;
%vector de ganancias
K = [k1 k2 k3 k4];
%referencia deseada 
ref = 4*pi;

%calculo del control obtenido del sistema 
v = -k1*(z1-ref)-k2*z2-k3*z3-k4*z4;
u = (I*J)/(k)*(v-a*sin(x1)*(x2^2+a*cos(x1)+b)-b*(x1-x3)*(b+c+a*cos(x1)));

%Sistema en la forma observador 
A = [0 1 0 0;-b 0 b 0;0 0 0 1;c 0 -c 0];

psi = [0;-a*sin(x1);0;d*u ];

C = [1 0 0 0];

xp = A*x+psi;

y = C*x;

%valores proporcionados
H = [149.1; 7440.6; 1391.5; 1714.6];
% H = [46;591;14.3; -419];
Xhp = A*xh +psi+H*(y-C*xh);

Xp = [xp;Xhp];

end
