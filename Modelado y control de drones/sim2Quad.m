clc 
clear all 
close all
%x y z phi theta psi xp yp zp phip thetap psip
T = 0.01;                       % periodo de simulación 
t = [0:T:20];                   % tiempo de simulación
x = [0 0 0 0 0 0 0 0 0 0 0 0];  % vector variables de estado                              
%----------------------condiciones iniciales-------------------------------
m = 1;
g = 9.81;
Ixx = 0.1;
Iyy = 0.1;                      
Izz = 0.1;
%--------------------------------ganancias---------------------------------
k1x = 4;
k2x = 12;
%------------------------------
k1y = 4;
k2y = 12;
%------------------------------
k1z = 6;
k2z = 25;
%------------------------------
k1phi = 7;
k2phi = 70; 
k1theta = 5; 
k2theta = 100;
k1psi = 7;
k2psi = 50;
k1 = [k1phi 0 0;0 k1theta 0;0 0 k1psi];
k2 = [k2phi 0 0;0 k2theta 0;0 0 k2psi];
%------------------------Dinámicas del error para x------------------------
xd = sin(t);            
xdp = cos(t);
xd2p = -sin(t);
%-----------------------Dinámicas del error para y-------------------------
yd = cos(t);            
ydp = -sin(t);
yd2p = -cos(t);
%---------------------- Dinámicas del error para z-------------------------
zd = 1;
zdp = 0;
zd2p = 0;
% zdp = -0.5*sin(t);
% zd2p = -0.5*cos(t);

%---------------------Dinámicas del error para psi-------------------------
psid = 0.5*cos(t);
%psidp = -0.5*sin(t);
%psid2p = -0.5*cos(t);
% psidp = 0;
% psid2p = 0;
%-------------------------entradas de control------------------------------
tau = [0 0 0]';
u = 0;
mux =0;                %Entrada de control virtual para x
muy =0;                %Entrada de control virtual para y

for i=1:length(t)-1            % length calcula la longitud del vector t
    [tt, xx] = ode45(@model2Quad,[t(i) t(i+1)],x(i,:),[],tau(:,i),u(i),mux(i),muy(i));
    x(i+1,:) = xx(end,:);
    %----------------Ecuaciones para la ley de control para fq-------------
    ez = x(i+1,3)-zd;
    ezp = x(i+1,9);
    u(i+1) = (m/(cos(x(i+1,4))*cos(x(i+1,5))))*(-k1z*ezp-k2z*ez+g);
    %-------------------phi theta psi phip thetap psip---------------------
    phi = x(i+1,4);
    theta = x(i+1,5);
    psi = x(i+1,6);
    phip = x(i+1,10);
    thetap = x(i+1,11);
    psip = x(i+1,12);
    etap =[phip thetap psip]';
    
    psid(i+1) = psid(i);
    %--------------Ecuaciones para la ley de control para mux--------------
    ex = x(i+1,1) - xd(i+1);
    e_xp = x(i+1,7) -xdp(i+1);
    mux(i+1) = (xd2p(i+1) - k1x*e_xp - k2x*ex);
    %-------------Ecuaciones para la ley de control para muy---------------
    ey = x(i+1,2) - yd(i+1);
    eyp = x(i+1,8) -ydp(i+1);
    muy(i+1) = (yd2p(i+1) - k1y*eyp - k2y*ey); 
    %------------Ecuaciones para el control deorientación------------------
    c11 = 0;
    c12 = (Iyy-Izz)*(thetap*cos(phi)*sin(phi)+psip*sin(phi)*sin(phi)...
        *cos(theta))+(Izz-Iyy)*psip*cos(phi)*cos(phi)*cos(theta)...
        -Ixx*psip*cos(theta);
    c13 = (Izz-Iyy)*psip*cos(phi)*sin(phi)*cos(theta)*cos(theta);
    c21 = (Izz-Iyy)*(thetap*cos(phi)*sin(phi)+psip*sin(phi)*sin(phi)...
        *cos(theta))+(Iyy-Izz)*psip*cos(phi)*cos(phi)*cos(theta)...
        +Ixx*psip*cos(theta);
    c22 = (Izz-Iyy)*phip*cos(phi)*sin(phi);
    c23 = -Ixx*psip*sin(theta)*cos(theta)+Iyy*psip*sin(phi)*sin(phi)...
        *sin(theta)*cos(theta)+Izz*psip*cos(phi)*cos(phi)*sin(theta)...
        *cos(theta);
    c31 = (Iyy-Izz)*psip*cos(theta)*cos(theta)*sin(phi)*cos(phi)...
        -Ixx*thetap*cos(theta);
    c32 = (Izz-Iyy)*(thetap*cos(phi)*sin(phi)*sin(theta)+phip*sin(phi)...
        *sin(phi)*cos(theta))+(Iyy-Izz)*thetap*cos(phi)*cos(phi)...
        *cos(theta)+Ixx*psip*sin(theta)*cos(theta)-Iyy*psip*sin(phi)...
        *sin(phi)*sin(theta)*cos(theta)-Izz*psip*cos(phi)*cos(phi)...
        *sin(theta)*cos(theta);
    c33 = (Iyy-Izz)*phip*cos(phi)*sin(phi)*cos(theta)*cos(theta)...
        -Iyy*thetap*sin(phi)*sin(phi)*cos(theta)*sin(theta)...
        -Izz*thetap*cos(phi)*cos(phi)*cos(theta)*sin(theta)...
        +Ixx*thetap*cos(theta)*sin(theta);
    C = [c11 c12 c13; c21 c22 c23; c31 c32 c33];
    J = [Ixx, 0, -Ixx*sin(theta); 
        0, Iyy*cos(phi)*cos(phi)+Izz*sin(phi)*sin(phi), ...
        cos(phi)*cos(theta)*sin(phi)*(Iyy-Izz); 
        -Ixx*sin(theta), cos(phi)*cos(theta)*sin(phi)*(Iyy-Izz), ...
        (Ixx)*sin(theta)*sin(theta)+Iyy*cos(theta)*cos(theta)*sin(phi)...
        *sin(phi)+Izz*cos(phi)*cos(phi)*cos(theta)*cos(theta)];
    %------------Ecuaciones para el seguimiento de trayectorias------------
    %-----de las dinámicas de X y Y es necesario que los ángulos-----------
    %----\phi y \theta sean los dados por las siguientes ecuaciones--------
    Phid(i+1)= asin((m/u(i+1))*(sin(x(i+1,6))*mux(i+1)-cos(x(i+1,6))*muy(i+1)));
    Thetad(i+1)= asin((m/u(i+1))*((cos(x(i+1,6))*mux(i+1)+sin(x(i+1,6))*muy(i+1)))/(cos(x(i+1,4))));
    %--------------Ecuaciones para la ley de control para Tau--------------
    error_eta = [Phid(i+1) Thetad(i+1) psid(i)]' - x(i+1,4:6)';
    error_etap = - x(i+1,10:12)'; 
    tau(:,i+1) = C*etap + J*( k2*error_eta + k1*error_etap);
end

%Graficas de orientación
figure(1)
subplot(311)
plot(t,Phid,t,x(:,4))
legend('\phi_d','\phi')
xlabel('Tiempo [seg]')
ylabel('Orientación en \phi [rad]')
subplot(312)
plot(t,Thetad,t,x(:,5))
legend('\theta_d','\theta')
xlabel('Tiempo [seg]')
ylabel('Orientación en \theta [rad]')
subplot(313)
plot(t,psid,t,x(:,6))
legend('\psi_d','\psi')
xlabel('Tiempo [seg]')
ylabel('Orientación en \psi [rad]')
%Graficas de posición
figure(2)
subplot(311)
plot(t,xd,t,x(:,1))
legend('X_d','X')
xlabel('Tiempo [seg]')
ylabel('Posición en x [m]')
subplot(312)
plot(t,yd,t,x(:,2))
legend('Y_d','Y')
xlabel('Tiempo [seg]')
ylabel('Posición en y [m]')
subplot(313)
plot(t,zd,t,x(:,3))
legend('Z_d','Z')
xlabel('Tiempo [seg]')
ylabel('Posición en z [m]')
%Gráfica 3D 
figure
plot3(x(:,1),x(:,2),x(:,3),xd,yd,1*ones(1,length(t)))
legend('Trayectoria del control','Trayectoria deseada')
grid minor
xlabel('Posición en x [m]')
ylabel('Posición en y [m]')
zlabel('Posición en z [m]')

