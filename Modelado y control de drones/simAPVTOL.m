clc
clear all
close all

T=0.01; %Definiendo los intervalos de tiempo
t = [0:T:10]; %vector de tiempos
x = [0 0 0.3 0 0 0];

%Entradas de control
fq = 0;
tauT = 0;
ux = 0; %Entrada de control virtual para x

%Dinámica del error para z
zd = 1 + 0.5*cos(t); 
zdp = -0.5*sin(t);
zd2p = -0.5*cos(t);
k1z = 6.1;
k2z = 7.5;

%Dinámica del error para x
xd = 1 + 0.5*cos(t);
xdp = -0.5*sin(t);
xd2p = -0.5*cos(t);
k1x = 4;
k2x = 4.5;

%Dinámica del error para theta
td = 0;
tdp = 0;
td2p = 0;
k1t = 6;
k2t = 120;

g = 9.81;
m = 1;
Iyy = 0.1;

for i=1:length(t)-1
    %1.- Medir las variables a controlar

    [tt,xx] = ode45(@modeloAPVTOL,[t(i) t(i+1)],x(i,:),[],ux(i),fq(i),tauT(i));
    x(i+1,:) = xx(end,:);

    %2.- Calcular la ley de control
    
    %Ley de control para z
    ez = x(i+1,2) - zd(i+1);
    ezp = x(i+1,5) - zdp(i+1);
    fq(i+1) = (m/cos(x(i+1,3)))*(zd2p(i+1)-k1z*ezp-k2z*ez+g);

    %Ley de control para x
    ex = x(i+1,1) - xd(i+1);
    exp = x(i+1,4) - xdp(i+1);
    ux(i+1) = m*(xd2p(i+1) - k1x*exp - k2x*ex);

    td(i+1) = asin(ux(i+1)/fq(i+1));
    
    %Ley de control para theta
    et = x(i+1,3) - td(i+1);
    etp = x(i+1,6) - tdp;
    tauT(i+1) = Iyy*(td2p-k1t*etp-k2t*et);

end

%Gráficas:

figure(1)
subplot(311)
plot(t,zd,t,x(:,2))
legend('z_d','z')
ylabel('Posición en Z[m]')

subplot(312)
plot(t,xd,t,x(:,1))
legend('x_d','x')
ylabel('Posición en X[m]')

subplot(313)
plot(t,td,t,x(:,3))
legend('\theta_d','\theta')
ylabel('Posición en \theta[rad]')
xlabel('Tiempo [seg]')

figure(2)
plot(t,fq)
ylabel('f_q [N]')
xlabel('Tiempo [s]')
