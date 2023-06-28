clc 
clear 
close all

T = 0.01;                        
t = [0:T:20];                   
deltaT = T;
x = [0 0 0 0 0 0 0 0 0 0 0 0];                                

m = 11.5;
W = 112.8;
B = 114.8;
rb = [0; 0; 0];
rg = [0; 0; 0.02];

xg = 0;
yg = 0;
zg = 0.02;

Ixy = 0;
Ixz = 0;
Iyz = 0;
Iyx = 0;
Izx = 0;
Izy = 0;

Ix = 0.16;
Iy = 0.16;
Iz = 0.16;

Xup = -5.5;
Yvp = -12.7;
Zwp = -14.57;
Kpp = -0.12;
Mqp = -0.12;
Nrp = -0.12;

Xu = -4.03;
Yv = -6.22;
Zw = -5.18;
Kp = -0.07;
Mq = -0.07;
Nr = -0.07;

XuNL = -18.18;
YvNL = -21.66;
ZwNL = -36.99;
KpNL = -1.55;
MqNL = -1.55;
NrNL = -1.55;

u = x(7);
v = x(8);
w = x(9);
p = x(10);
q = x(11);
r = x(12);

eta = [x(1); x(2); x(3); x(4); x(5); x(6)];
nu = [x(7); x(8); x(9); x(10); x(11); x(12)];

phi = x(4);
theta = x(5);
psi = x(6);

kp = [3 3 3 4 4 2];
ki = [0.2 0.2 0.2 0.3 0.3 0.1];
kd = [2.5 2.5 0.5 0.5 1 0.5];

TauPID = [0 0 0 0 0 0]';
% etad = [xd; yd; zd; phid; thetad; psid];

for i=1:length(t)-1            
    [tt, xx] = ode45(@modelROV,[t(i) t(i+1)],x(i,:),[],TauPID(:,i));
    x(i+1,:) = xx(end,:);

    J = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(theta)*sin(psi), 0,           0,             0;
        cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi),  cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi),   0,           0,             0;
        -sin(theta),            sin(phi)*cos(theta),                                        cos(phi)*cos(theta),                0,           0,             0;
            0,                              0,                                                      0,                          1, sin(phi)*tan(theta), cos(phi)*tan(theta);
            0,                              0,                                                      0,                          0,       cos(phi),          -sin(phi);
            0,                              0,                                                      0,                          0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
    M_RB = [m    0    0    0    m*zg -m*yg;
            0    m    0  -m*zg   0    m*xg;
            0    0    m   m*yg -m*xg   0;
            0  -m*zg m*yg  Ix  -Ixy  -Ixz;
            m*zg  0   -m*xg -Iyx  Iy   -Iyz;
            -m*yg  m*xg  0   -Izx -Izy   Iz];
    M_A = -[Xup 0   0   0   0   0;
            0   Yvp 0   0   0   0;
            0   0  Zwp  0   0   0;
            0   0   0   Kpp 0   0;
            0   0   0   0   Mqp 0;
            0   0   0   0   0   Nrp];    
    C_RB = [0    0    0     0    m*w   0;
            0    0    0   -m*w   0     0;
            0    0    0    m*v  -m*u   0;
            0    m*w -m*v   0    Iz*r -Iy*q;
          -m*w  0   -m*u -Iz*r   0    Ix*p;
           m*v -m*u  0    Iy*q -Ix*p  0];
    C_A = [0     0      0       0   Zwp*w  0;
           0     0      0    -Zwp*w  0    -Xup*u;
           0     0      0    -Yvp*v Xup*u  0;
           0    -Zwp*w Yvp*v   0   -Nrp*r -Mqp*q;
          Zwp*w 0    -Xup*u  Nrp*r  0    -Kpp*p
         -Yvp*v Xup*u  0    -Mqp*q Kpp*p  0];
    M = M_RB + M_A;
    C = C_RB + C_A;
    D = [Xu+XuNL    0       0   0       0       0;
            0   Yv+YvNL     0   0       0       0;
            0       0   Zw+ZwNL 0       0       0;
            0       0       0   Kp+KpNL 0       0;
            0       0       0   0       Mq+MqNL 0;
            0       0       0   0       0       Nr+NrNL];
    g = [(W-B)*sin(theta); -(W-B)*cos(theta)*sin(phi); -(W-B)*cos(theta)*cos(phi); zg*W*cos(theta)*sin(phi); zg*W*sin(theta);0];
    
    etad = [0; 0; 10; 0; 0; 0];

    xd(i+1) = etad(1);
    yd(i+1) = etad(2);
    zd(i+1) = etad(3);
    phid(i+1) = etad(4);
    thetad(i+1) = etad(5);
    psid(i+1) = etad(6);
    
    etap = J*nu;
    ep = -(J*nu);

    xsim =x(i+1,1);
    ysim =x(i+1,2);
    zsim =x(i+1,3);
    phi2 =x(i+1,4);
    theta2 =x(i+1,5);
    psi2 =x(i+1,6);
    
    eta = [xsim ysim zsim phi2 theta2 psi2];
    
    e = etad - eta;
    eb = J.'*e.';
    epb = J.'*ep;
    eib =  eb + T*epb;
    TauPID(:,i+1) = kp*eb + ki*eib + kd*epb;
    
    if etad == [0; 0; 10; 0; 0; 0]
        [tt, xx] = ode45(@modelROV,[t(i) t(i+1)],x(i,:),[],TauPID(:,i));
        x(i+1,:) = xx(end,:);
    elseif etad == [0; 5; 10; 0; 0; 0]
        [tt, xx] = ode45(@modelROV,[t(i) t(i+1)],x(i,:),[],TauPID(:,i)); 
        x(i+1,:) = xx(end,:);
    elseif etad == [5; 5; 10; 0; 0; 0]
        [tt, xx] = ode45(@modelROV,[t(i) t(i+1)],x(i,:),[],TauPID(:,i));
        x(i+1,:) = xx(end,:);
    elseif etad == [5; 0; 10; 0; 0; 0]
        [tt, xx] = ode45(@modelROV,[t(i) t(i+1)],x(i,:),[],TauPID(:,i));
        x(i+1,:) = xx(end,:);
    elseif etad == [0; 0; 10; 0; 0; 0]
        [tt, xx] = ode45(@modelROV,[t(i) t(i+1)],x(i,:),[],TauPID(:,i));
        x(i+1,:) = xx(end,:);
    end
    
end

figure(1)
plot(t,xd, t, x(:,1))
legend('trayectoria deseada', 'trayectoria simulada')
xlabel('Tiempo [seg]')
ylabel('posición x [m]')

figure(2)
plot(t,yd, t, x(:,2))
legend('trayectoria deseada', 'trayectoria simulada')
xlabel('Tiempo [seg]')
ylabel('posición y [m]')

figure(3)
plot(t,zd, t, x(:,3))
legend('trayectoria deseada', 'trayectoria simulada')
xlabel('Tiempo [seg]')
ylabel('posición z [m]')

figure(4)
plot(t,phid, t, x(:,4))
legend('orientación deseada', 'orientación simulada')
ylabel('orientación del eje x [rad]')
xlabel('Tiempo [seg]')

figure(5)
plot(t,thetad, t, x(:,5))
legend('orientación deseada', 'orientación simulada')
ylabel('orientación del eje y [rad]')
xlabel('Tiempo [seg]')

figure(6)
plot(t,psid, t, x(:,6))
legend('orientación deseada', 'orientación simulada')
ylabel('orientación del eje z [rad]')
xlabel('Tiempo [sef]')

figure(7)
plot(t,xd, t, x(:,1))
plot(t,zd, t, x(:,3))
hold on
xlabel('posición x [m]')
ylabel('posición z [m]')

figure (8)
plot(x(:,1),x(:,2),[0, 5], [0,0], 'g',[5, 5], [0,5], 'g', [0,5], [5, 5], 'g', [0 ,0], [5,0], 'g')
legend('Trayectoria simulada','Trayectoria deseada')
xlabel('posición x [m]')
ylabel('posición y [m]')

figure(9)
plot(t, x(:,1), t, x(:,3))
xlabel('posición x [m]')
ylabel('posición y [m]')

%Gráfica 3D
figure(10)
plot3(zd,yd,2*ones(1,length(t)),x(:,3),x(:,2),2*ones(1,length(t)))
legend('Trayectoria deseada','Trayectoria simulada')
xlabel('posición x [m]')
ylabel('posición y [m]')
zlabel('posición z [m]')
grid minor 

figure(11)
plot3(x(:,1),x(:,2),x(:,3),xd,yd,2*ones(1,length(t)))
legend('Trayectoria simulada','Trayectoria deseada')
grid minor
xlabel('Posición en x [m]')
ylabel('Posición en y [m]')
zlabel('Posición en z [m]')


