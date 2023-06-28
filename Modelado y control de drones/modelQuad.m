function dxdt = modelQuad(t,x,tau,f)

m = 1;
g = 9.81;
Ixx = 0.1;
Iyy = 0.1;
Izz = 0.1;
phi = x(4);
theta = x(5);
psi = x(6);
phip = x(10);
thetap = x(11);
psip = x(12);

etap = [phip;thetap;psip];

c11 = 0;
c12 = (Iyy-Izz)*(thetap*cos(phi)*sin(phi)+psip*sin(phi)*sin(phi)*cos(theta))+(Izz-Iyy)*psip*cos(phi)*cos(phi)*cos(theta)-Ixx*psip*cos(theta);
c13 = (Izz-Iyy)*psip*cos(phi)*sin(phi)*cos(theta)*cos(theta);
c21 = (Izz-Iyy)*(thetap*cos(phi)*sin(phi)+psip*sin(phi)*sin(phi)*cos(theta))+(Iyy-Izz)*psip*cos(phi)*cos(phi)*cos(theta)+Ixx*psip*cos(theta);
c22 = (Izz-Iyy)*phip*cos(phi)*sin(phi);
c23 = -Ixx*psip*sin(theta)*cos(theta)+Iyy*psip*sin(phi)*sin(phi)*sin(theta)*cos(theta)+Izz*psip*cos(phi)*cos(phi)*sin(theta)*cos(theta);
c31 = (Iyy-Izz)*psip*cos(theta)*cos(theta)*sin(phi)*cos(phi)-Ixx*thetap*cos(theta);
c32 = (Izz-Iyy)*(thetap*cos(phi)*sin(phi)*sin(theta)+phip*sin(phi)*sin(phi)*cos(theta))+(Iyy-Izz)*thetap*cos(phi)*cos(phi)*cos(theta)+Ixx*psip*sin(theta)*cos(theta)-Iyy*psip*sin(phi)*sin(phi)*sin(theta)*cos(theta)-Izz*psip*cos(phi)*cos(phi)*sin(theta)*cos(theta);
c33 = (Iyy-Izz)*phip*cos(phi)*sin(phi)*cos(theta)*cos(theta)-Iyy*thetap*sin(phi)*sin(phi)*cos(theta)*sin(theta)-Izz*thetap*cos(phi)*cos(phi)*cos(theta)*sin(theta)+Ixx*thetap*cos(theta)*sin(theta);
C = [c11 c12 c13; c21 c22 c23; c31 c32 c33];
J = [Ixx, 0, -Ixx*sin(theta); 0, Iyy*cos(phi)*cos(phi)+Izz*sin(phi)*sin(phi), cos(phi)*cos(theta)*sin(phi)*(Iyy-Izz); -Ixx*sin(theta), cos(phi)*cos(theta)*sin(phi)*(Iyy-Izz), (Ixx)*sin(theta)*sin(theta)+Iyy*cos(theta)*cos(theta)*sin(phi)*sin(phi)+Izz*cos(phi)*cos(phi)*cos(theta)*cos(theta)];


%------------------Ecuaciones representación de estado-----------------------
%---------------x y z phi theta psi xp yp zp phip thetap psip----------------

dxdt(1,1) = x(7);
dxdt(2,1) = x(8);
dxdt(3,1) = x(9);
dxdt(4,1) = x(10);
dxdt(5,1) = x(11);
dxdt(6,1) = x(12);
dxdt(7,1) = (f/m)*(cos(x(4))*sin(x(5))*cos(x(6))+sin(x(4))*sin(x(6)));
dxdt(8,1) = (f/m)*(cos(x(4))*sin(x(5))*sin(x(6))-sin(x(4))*cos(x(6)));
dxdt(9,1) = (f/m)*cos(x(4))*cos(x(5))-g;
dxdt(10:12,1) = inv(J)*(tau - C*etap);







