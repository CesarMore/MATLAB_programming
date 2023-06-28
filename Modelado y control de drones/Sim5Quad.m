clc
clear all
close all

%%%En este programa se controla solamente x y z phi theta y psi %%%

T = 0.01;
t = [0:T:10];
x = [0 0 0 0 0 0 0 0 0 0 0 0];
u = 0;
Tau = [0;0;0];
Mux = 0;
Muy = 0;
PHID = 0;
PRUEBAPHID = 0;

g = 9.81;
m = 1;
Ixx = 0.1;
Iyy = 0.1;
Izz = 0.1;

zd = 1;
zdp = 0;
zd2p = 0;
xd = cos(t);
xdp = -sin(t);
xd2p = -cos(t);
yd = sin(t);
ydp = cos(t);
yd2p = -sin(t);

% phid = +0.5*cos(t);
% phidp = -0.5*sin(t);
% phid2p = -0.5*cos(t);
% td = +0.5*cos(t);
% tdp = -0.5*sin(t);
% td2p = -0.5*cos(t);
% psid = +0.5*cos(t);
psid = 0;
% psidp = -0.5*sin(t);
% psid2p = -0.5*cos(t);

% k1z = 120.1;
% k2z = 350.5;
k1z = 10.1;
k2z = 30.5;
k1x = 6;
k2x = 60;
k1y = 6;
k2y = 40;


% k1phi = 10; 
% k1th = 10.5;
% k1psi = 3.5;
% k2phi = 120; 
% k2th = 120;
% k2psi = 20;
k1phi = 7;
k2phi = 50; 
k1th = 7; 
k2th = 60;
k1psi = 7;
k2psi = 50;
k1T = [k1phi 0 0;0 k1th 0;0 0 k1psi];
k2T = [k2phi 0 0;0 k2th 0;0 0 k2psi];

for i=1:length(t)-1
    [tt,xx] = ode45(@modelo3Quad,[t(i) t(i+1)],x(i,:),[],Tau(:,i),Mux(i),Muy(i),u(i));
    x(i+1,:) = xx(end,:);
     
    %zd(i+1) = zd(i);
    ez = x(i+1,3)-zd;
    ezp = x(i+1,9);
    u(i+1) = (m/(cos(x(i+1,4))*cos(x(i+1,5))))+(-k1z*ezp-k2z*ez+g);

    phi = x(i+1,4);
    theta = x(i+1,5);
    psi = x(i+1,6);
    phip = x(i+1,10);
    thetap = x(i+1,11);
    psip = x(i+1,12);
    etap = [phip;thetap;psip];
    psid(i+1) = psid(i);
    
    ex = x(i+1,1)-xd(i+1);
    exp = x(i+1,7)-xdp(i+1);
    Mux(i+1) = (xd2p(i+1)-k1x*exp-k2x*ex);
    
    ey = x(i+1,2)-yd(i+1);
    eyp = x(i+1,8)-ydp(i+1);
    Muy(i+1) = (yd2p(i+1)-k1y*eyp-k2y*ey);
    
    Jota = [Ixx 0 -Ixx*sin(theta);
        0 Iyy*cos(phi)^2+Izz*sin(phi)^2 cos(phi)*cos(theta)*sin(phi)*(Iyy-Izz);
        -Ixx*sin(theta) cos(phi)*cos(theta)*sin(phi)*(Iyy-Izz) Ixx*sin(theta)^2+Iyy*(cos(theta)^2)*(sin(phi)^2)+Izz*(cos(phi)^2)*(cos(theta)^2)];

    C = [0 (Iyy-Izz)*(thetap*cos(phi)*sin(phi)+psip*sin(phi)^2*cos(theta))+(Izz-Iyy)*(psip*cos(phi)^2*cos(theta))-Ixx*psip*cos(theta) (Izz-Iyy)*psip*cos(phi)*sin(phi)*cos(theta)^2;
        (Izz-Iyy)*(thetap*cos(phi)*sin(phi)+psip*sin(phi)^2*cos(theta))+(Iyy-Izz)*psip*cos(phi)^2*cos(theta)+Ixx*psip*cos(theta) (Izz-Iyy)*phip*cos(phi)*cos(theta) -Ixx*psip*sin(theta)*cos(theta)+Iyy*psip*sin(phi)^2*sin(theta)*cos(theta)+Izz*psip*cos(phi)^2*sin(theta)*cos(theta);
        (Iyy-Izz)*(psip*sin(phi)*cos(phi)*cos(theta)^2)-Ixx*thetap*cos(theta) (Izz-Iyy)*(thetap*sin(phi)*cos(phi)*sin(theta)+phip*sin(phi)^2*cos(theta))+(Iyy-Izz)*(phip*cos(phi)^2*cos(theta))+Ixx*psip*sin(theta)*cos(theta)-Iyy*psip*sin(phi)^2*sin(theta)*cos(theta)-Izz*psip*cos(phi)^2+sin(theta)*cos(theta) (Iyy-Izz)*(phip*sin(phi)*cos(phi)*cos(theta)^2)-Iyy*thetap*sin(phi)^2*sin(theta)*cos(theta)-Izz*thetap*cos(phi)^2*sin(theta)*cos(theta)+Ixx*thetap*sin(theta)*cos(theta)];
    
    
    
    %PRUEBAPHID(i+1) = (m/u(i+1))*(sin(psid(i+1))*Mux(i+1)-cos(psid(i+1))*Muy(i+1))
    
    %%%Tome psi y phi calculadas, no deseadas
    PHID(i+1) = asin((m/u(i+1))*(sin(x(i+1,6))*Mux(i+1)-cos(x(i+1,6))*Muy(i+1)));
    
    THETAD(i+1) = asin((m/u(i+1))*((cos(x(i+1,6))*Mux(i+1)+sin(x(i+1,6))*Muy(i+1)))/(cos(x(i+1,4))));
    
    epunto = [0-x(i+1,10);0-x(i+1,11);0-x(i+1,12)];
    e = [PHID(i+1)-x(i+1,4);THETAD(i+1)-x(i+1,5);psid(i+1)-x(i+1,6)];
    
    %eta2pd = [phid2p(i+1);td2p(i+1);psid2p(i+1)];
    
    Tau(:,i+1) = C*etap+Jota*(k1T*epunto+k2T*e);
    
end


figure(1)
plot(t,zd,t,x(:,3))
legend('zd','z')

figure(2)
plot(t,PHID,t,x(:,4))
legend('phid','phi')

figure(3)
plot(t,THETAD,t,x(:,5))
legend('td','t')

figure(4)
plot(t,psid,t,x(:,6))
legend('psid','psi')

figure(5)
plot(t,Tau(1,:))
legend('Tauphi')
figure(6)
plot(t,Tau(2,:))
legend('Tautheta')
figure(7)
plot(t,Tau(3,:))
legend('Taupsi')

figure(8)
subplot(211)
plot(t,Mux)
legend('Mux')
subplot(212)
plot(t,Muy)
legend('Muy')

figure(9)
plot(t,PHID)
figure(10)
plot(t,PRUEBAPHID)

figure(11)
plot(t,u)
legend('u')

figure(12)
plot(t,xd,t,x(:,1))
legend('xd','x')

figure(13)
plot(t,yd,t,x(:,2))
legend('yd','y')

figure(14)
plot3(x(:,1),x(:,2),x(:,3))
grid on

figure(15)
plot3(x(:,1),x(:,2),x(:,3),xd,yd,1*ones(1,length(t)))
legend('Trayectoria de control','Trayectoria deseada')
grid minor

% figure(9)
% plot(t,xd,t,x(:,1))
% legend('xd','x')
% 
% figure(10)
% plot(t,yd,t,x(:,2))
% legend('yd','y')
