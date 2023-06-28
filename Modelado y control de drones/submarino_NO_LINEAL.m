clc 
clear 
close all
%-------------Tiempos para cada trayectoria-------------------                       
T = 0.05;                      %Intervalo de muestreo 
t = [0:T:80];     % Tiempo de la trayectoria de ascenso  

t2= [0:T:80];     %Tiempo de la trayectoria del cuadrado
%-------------------------------------------------------------------------
%-------------------------- DATOS EXPERIMENTALES -------------------------
W=112.8;    %m*g (N)
B=114.8;    % (N)
m=11.5;     %(kg)
%rb=[0;0;0];     Vector de posición del centro geometrico
%rg=[0;0;0.02]   Vector de posición del centro de gravedad
xg=0;
yg=0;
zg=0.02;
Ix = 0.16;  %(kg m^2)
Iy = 0.16;  %(kg m^2)
Iz = 0.16;  %(kg m^2)
%----------------------
Xup=-5.5;        %(kg)
Yvp=-12.7;       %(kg)
Zwp=-14.57;      %(kg)
Kppunto=-0.12;   %(kg m^2/rad)
Mqp=-0.12;       %(kg m^2/rad)
Nrp=-0.12;       %(kg m^2/rad)
%------------------------
Xu=-4.03;    %(Ns/m)
Yv=-6.22;    %(Ns/m)
Zw=-5.18;    %(Ns/m)
Kp=-0.07;    %(Ns/rad)
Mq=-0.07;    %(Ns/rad)
Nr=-0.07;    %(Ns/rad)
%------------------------
Xuu=-18.18;    %(Ns^2/m^2)
Yvv=-21.66;    %(Ns^2/m^2)
Zww=-36.99;    %(Ns^2/m^2)
Kpp=-1.55;     %(Ns^2/rad^2)
Mqq=-1.55;     %(Ns^2/rad^2)
Nrr=-1.55;     %(Ns^2/rad^2)
%-------------------------------------------------
%--------------  Valores iniciales de x para cada trayecto -------------
x = [0 0 0 0 0 0 0 0 0 0 0 0];         % V.I.  ascenso 
%x2= [0 0 x(end,3) 0 0 0 0 0 0 0 0 0]; %  V.I. cuadrado

%------------------- Eta deseada para cada trayectoria --------------------
%---------------------- etad=[x  y  z  phi theta psi] ---------------------
etad=[ 0 0 10 0 0 0]'; %eta deseada ascenso
etad2=[5 0 10 0 0 0]'; %eta deseada cuadrado

%-------------------------- TauPID inicial -------------------------------
% -------------------- tauPID=[ X ;Y ;Z ;K ;M ;N]-------------------------
tauPID = [0; 0; 0; 0; 0; 0];  % Trayectoria ascenso
tauPID2= [0; 0; 0; 0; 0; 0];  % Trayectoria cuadrado                     
%------------------------ TAU NO LINEAL-----------------------------------
tauPIDNL = [0; 0; 0; 0; 0; 0];  % Trayectoria ascenso NO LINEAL
tauPID2NL= [0; 0; 0; 0; 0; 0];  % Trayectoria cuadrado  NO LINEAL
%--------------------------------Ganancias---------------------------------
%------------------------------ PID NO LINEAL -----------------------------
KP=diag([16           8        28      30       10        50]);
KI=diag([11           25       10       1        5       100]);
KD=diag([41           26       18       0        1       55 ]);
              
%--------------------------------------------------------------------------
%----------------------- Trayectoria de ascenso ---------------------------
%--------------------------------------------------------------------------
for i=1:length(t)-1            
    [tt, xx] = ode45(@modeloSubNL,[t(i) t(i+1)],x(i,:),[],tauPIDNL(:,i));
    x(i+1,:) = xx(end,:);
%--------------------------------------------------------------------------
    xreal=x(i+1,1);
    yreal=x(i+1,2);
    zreal=x(i+1,3);
    phi = x(i+1,4);
    theta = x(i+1,5);
    psi = x(i+1,6);
%------------------------
     nu= x(i+1,7:12)';
% ---------------- Cambio de marco de referencia  MATRIZ J ---------------
J = [cos(psi)*cos(theta)    -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi)      sin(psi)*sin(phi)+cos(psi)*sin(theta)*cos(phi)     0           0                   0;  
     sin(psi)*cos(theta)    cos(psi)*cos(phi)+sin(psi)*sin(theta)*sin(phi)        -cos(psi)*sin(phi)+sin(psi)*sin(theta)*cos(phi)   0           0                   0;
     -sin(theta)                         cos(theta)*sin(phi)                                 cos(theta)*cos(phi)                    0           0                   0;
         0                                      0                                                       0                           1   tan(theta)*sin(phi)   tan(theta)*cos(phi);
         0                                      0                                                       0                           0   cos(phi)              -sin(phi);
         0                                      0                                                       0                           0   sin(phi)/cos(theta)   cos(phi)/cos(theta)]; 
JT= transpose(J); % J TRAMPUESTA 

%---------------------------- error eta  ----------------------------------
%error eta = eta deseada - eta
% etad(i+1,:)= [xd yd zd phid thetad psid]';
eeta= etad -x(i+1,1:6)';                   % INERCIAL 
eetab(:,i+1)=JT*eeta;                      % CUERPO   

%-------------------------- error eta punto -------------------------------
%error eta punto  = eta deseada punto - eta punto
%etap=J*nu (ecuacion 4.1) NU SALE DEL ODE45
etap=J*x(i+1,7:12)';
eetap=  - etap;                             % INERCIAL
eetapb(:,i+1)=JT*eetap;                     % CUERPO                                   
%-------------------------Integral del error ------------------------------
% Integración por Euler
% delta = T
     ebi(:,i+1)=eetab(:,i)+ T*eetapb(:,i);      
%--------------------------- PID LINEAL -----------------------------------
%-------------------- Ley de control tau LINEAL ---------------------------
    
tauPID(:,i+1)=KP*eetab(:,i+1)+KI*ebi(:,i+1)+KD*eetapb(:,i+1);
%--------------------------------------------------------------------------
phi = x(4);
theta = x(5);
psi = x(6);
u= x(7);
v= x(8);
w= x(9);
p= x(10);
q= x(11);
r= x(12);

%-------------------------- PID  NO LINEAL --------------------------------
MRB=[ m      0       0      0      m*zg   0;
      0      m       0      -m*zg  0      0;
      0      0       m      0      0      0;
      0      -m*zg   0      Ix     0      0;
      m*zg   0       0      0      Iy     0;
      0      0       0      0      0      Iz];
  
 
CRB=[0      0       0       0       m*w     0;
     0      0       0       -m*w    0       0;
     0      0       0       m*v     -m*u    0;
     0      m*w     -m*v    0       Iz*r    -Iy*q;
     -m*w   0       -m*u    -Iz*r   0       Ix*p;
     m*v    -m*u    0       Iy*q    -Ix*p   0 ];
 
MA=-[Xup 0  0  0  0  0;
     0  Yvp 0  0  0  0;
     0   0 Zwp 0  0  0;
     0   0  0 Kpp 0  0;
     0   0  0  0 Mqp 0;
     0   0  0  0  0  Nrp];

CA=[0          0       0        0         Zwp*w    0;
    0          0       0        -Zwp*w    0        -Xup*u;
    0          0       0        -Yvp*v    Xup*u    0;
    0          -Zwp*w  Yvp*v    0         -Nrp*r   Mqp*q;
    Zwp*w      0       -Xup*u   Nrp*r     0        -Kppunto*p;
    -Yvp*v     Xup*u   0        -Mqp*q    Kppunto*p    0];

D=-[Xu+Xuu*abs(u)   0               0                  0                 0                 0;
    0                Yv+Yvv*abs(v)  0                  0                 0                 0;
    0                0                Zw+Zww*abs(w)    0                 0                 0;
    0                0                0                  Kp+Kpp*abs(p)   0                 0;
    0                0                0                  0                 Mq+Mqq*abs(q)   0;
    0                0                0                  0                 0                 Nr+Nrr*abs(r)];

g=[(W-B)*sin(theta);
    -(W-B)*cos(theta)*sin(phi);
    -(W-B)*cos(theta)*cos(phi);
    zg*W*cos(theta)*sin(phi);
    zg*W*sin(theta);
    0];
%------------------------------------------------
M=MRB+MA;
C=CRB+CA;
%------------------- Ley de control tau NO LINEAL -------------------------
tauPIDD=tauPID(:,end);
tauPIDNL(:,i+1)= M*tauPIDD+C*nu+D*nu+g;
end

%--------------------------------------------------------------------------
%---------------------- Trayectoria del cuadrado --------------------------
%--------------------------------------------------------------------------
x2= [0 0 x(end,3) 0 0 0 0 0 0 0 0 0]; %Nuevos V.I.

for i=1:length(t2)-1            
    [tt2, xx2] = ode45(@modeloSubNL,[t2(i) t2(i+1)],x2(i,:),[],tauPID2NL(:,i));
    x2(i+1,:) = xx2(end,:);
%--------------------------------------------------------------------------
    xreal2=x2(i+1,1);
    yreal2=x2(i+1,2);
    zreal2=x2(i+1,3);
    phi2 = x2(i+1,4);
    theta2 = x2(i+1,5);
    psi2 = x2(i+1,6);
%------------------------
    nu2= x2(i+1,7:12)';
% ---------------- Cambio de marco de referencia  MATRIZ J ---------------
J2 = [cos(psi2)*cos(theta2)    -sin(psi2)*cos(phi2)+cos(psi2)*sin(theta2)*sin(phi2)      sin(psi2)*sin(phi2)+cos(psi2)*sin(theta2)*cos(phi2)     0           0                   0;  
     sin(psi2)*cos(theta2)    cos(psi2)*cos(phi2)+sin(psi2)*sin(theta2)*sin(phi2)        -cos(psi2)*sin(phi2)+sin(psi2)*sin(theta2)*cos(phi2)   0           0                   0;
     -sin(theta2)                         cos(theta2)*sin(phi2)                                 cos(theta2)*cos(phi2)                    0           0                   0;
         0                                      0                                                       0                           1   tan(theta2)*sin(phi2)   tan(theta2)*cos(phi2);
         0                                      0                                                       0                           0   cos(phi2)              -sin(phi2);
         0                                      0                                                       0                           0   sin(phi2)/cos(theta2)   cos(phi2)/cos(theta2)]; 
JT2= transpose(J2); % J2 TRAMPUESTA 

%---------------------------- error eta  ----------------------------------
%error eta = eta deseada - eta
% etad(i+1,:)= [xd yd zd phid thetad psid]';

% eeta= etad -x(i+1,1:6)';                    % INERCIAL 
% eetab(:,i+1)=JT*eeta;                       % CUERPO
eeta2= etad2 -x2(i+1,1:6)';                   % INERCIAL 
eetab2(:,i+1)=JT2*eeta2;                      % CUERPO

%-------------------------- error eta punto -------------------------------
%error eta punto  = eta deseada punto - eta punto
%etap=J*nu (ecuacion 4.1) NU SALE DEL ODE45
etap2=J2*x2(i+1,7:12)';
eetap2=  - etap2;                              % INERCIAL
eetapb2(:,i+1)=JT2*eetap2;                     % CUERPO                                   
%-------------------------Integral del error ------------------------------
% Integración por Euler
% delta = T
     ebi2(:,i+1)=eetab2(:,i)+ T*eetapb2(:,i);      
%--------------------------- PID LINEAL -----------------------------------
%--------------------- Ley de control tau ---------------------------------
    
tauPID2(:,i+1)=KP*eetab2(:,i+1)+KI*ebi2(:,i+1)+KD*eetapb2(:,i+1);
%--------------------------------------------------------------------------
phi2 = x2(4);
theta2 = x2(5);
psi2 = x2(6);
u2= x2(7);
v2= x2(8);
w2= x2(9);
p2= x2(10);
q2= x2(11);
r2= x2(12);

%-------------------------- PID  NO LINEAL --------------------------------
MRB2=[ m      0       0      0      m*zg   0;
      0      m       0      -m*zg  0      0;
      0      0       m      0      0      0;
      0      -m*zg   0      Ix     0      0;
      m*zg   0       0      0      Iy     0;
      0      0       0      0      0      Iz];
  
 
CRB2=[0      0        0        0        m*w2     0;
     0       0        0        -m*w2    0        0;
     0       0        0        m*v2     -m*u2    0;
     0       m*w2     -m*v2    0        Iz*r2    -Iy*q2;
     -m*w2   0        -m*u2    -Iz*r2   0        Ix*p2;
     m*v2    -m*u2    0        Iy*q2    -Ix*p2   0 ];
 
MA2=-[Xup 0  0  0  0  0;
     0  Yvp 0  0  0  0;
     0   0 Zwp 0  0  0;
     0   0  0 Kpp 0  0;
     0   0  0  0 Mqp 0;
     0   0  0  0  0  Nrp];

CA2=[0         0        0         0          Zwp*w2    0;
    0          0        0         -Zwp*w2    0         -Xup*u2;
    0          0        0         -Yvp*v2    Xup*u2    0;
    0          -Zwp*w2  Yvp*v2    0          -Nrp*r2   Mqp*q2;
    Zwp*w2     0        -Xup*u2   Nrp*r2     0         -Kpp*p2;
    -Yvp*v2    Xup*u2   0         -Mqp*q2    Kpp*p2    0];

D2=-[Xu+Xuu*abs(u2)   0               0                  0                 0                 0;
    0                Yv+Yvv*abs(v2)  0                  0                 0                 0;
    0                0                Zw+Zww*abs(w2)    0                 0                 0;
    0                0                0                  Kp+Kpp*abs(p2)   0                 0;
    0                0                0                  0                 Mq+Mqq*abs(q2)   0;
    0                0                0                  0                 0                 Nr+Nrr*abs(r2)];

g2=[(W-B)*sin(theta2);
    -(W-B)*cos(theta2)*sin(phi2);
    -(W-B)*cos(theta2)*cos(phi2);
    zg*W*cos(theta2)*sin(phi2);
    zg*W*sin(theta2);
    0];
%------------------------------------------------
M2=MRB2+MA2;
C2=CRB2+CA2;
%------------------- Ley de control tau NO LINEAL -------------------------
tauPID22=tauPID2(:,end);
tauPID2NL(:,i+1)= M2*tauPID22+C2*nu2+D2*nu2+g2;
%--------------------------------------------------------------------------
     if (xreal2 >= 5 )                % Para ir al punto 3
        etad2 = [5 5 10 0 0 0]';
     end  
     if (yreal2 >= 5)                 % Para ir al punto 4
        etad2 = [0 5 10 0 0 0]';
     end   
     if (xreal2 < 0)                  % Para ir al punto 1
        etad2 = [0 0 10 0 0 0]';
     end
end  
%--------------------------

%Grafica de trayectoria cuadrado
%Grafica x,y en 2D
figure(2)
plot(x2(:,1),x2(:,2),[0, 5], [0, 0],'g',[5, 5], [0, 5],'g',[0, 5], [5, 5],'g',[0, 0], [5, 0],'g')
legend('Trayectoria real','Trayectoria deseada')
ylabel('Posición en Y[m]')
xlabel('Posición en X[m]')

%Grafica trayectoria completa
% Grafica de trayectoria 3D
figure(3)
plot3(x2(:,1),x2(:,2),x2(:,3),x(:,1),x(:,2),x(:,3),[0,0],[0,0],[0,10],'g--',[0, 5], [0, 0],[10, 10],'g',[5, 5], [0, 5],[10, 10],'g',[0, 5], [5, 5],[10, 10],'g',[0, 0], [5, 0],[10, 10],'g')
xlim([-1 6])
ylim([-1 6])
zlim([0 11])
legend('Trayectoria para el  cuadrado','Trayectoria de ascenso','Trayectoria deseada')
xlabel('Posición en X[m]')
ylabel('Posición en Y[m]')
zlabel('Posición en Z[m]')
grid minor

%Grafica de posición cuadrado
figure(4)
subplot(311)
plot(t,x(:,1))
legend('X')
ylabel('Posición en X[m]')
subplot(312)
plot(t,x(:,2))
legend('Y')
ylabel('Posición en Y[m]')
subplot(313)
plot(t,x(:,3))
legend('Z')
ylabel('Posición en Z[m]')
xlabel('Tiempo[segundos]')
%-------------------------
%Grafica de posición ascenso
figure(5)
subplot(311)
plot(t,x2(:,1))
legend('X')
ylabel('Posición en X[m]')
subplot(312)
plot(t,x2(:,2))
legend('Y')
ylabel('Posición en Y[m]')
subplot(313)
plot(t,x2(:,3))
legend('Z')
ylabel('Posición en Z[m]')
xlabel('Tiempo[segundos]')
%-------------------------