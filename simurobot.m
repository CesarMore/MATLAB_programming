clear all
close all
clc

ti = 0;   %tiempo de integración
h = 0.01; %paso de integración 
tf = 5;   % tiempo final de integración 

ts = (ti:h:tf)'; %tiempo de simulación

x0 = [pi;0;0;0]; %valores inciales de los estados 
xh0 = [-0.8;0.1;0.1;0.1]; %valores iniciales de los estados calculados

config = odeset('RelTol',h,'AbsTol',h,'InitialStep',h,'MaxStep',h);
[t,x] = ode45('robotmani',ts,[x0;xh0],config);

%Estados calculados
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
%Estados estimados
xh1 = x(:,5);
xh2 = x(:,6);
xh3 = x(:,7);
xh4 = x(:,8);

% grafica del control u
[n,m] = size(t);
uu = zeros(n,1);
for k =1:n
    [Xp,u] = robot(t(k),[x1(k);x2(k);x3(k);x4(k);xh1(k);xh2(k);xh3(k);xh4(k)]);
    uu(k) = u;


end



%Gráfica de los estados vs estado estimado
subplot(2,2,1), plot(t,x1,t,xh1), legend('trayectoria deseada', 'trayectoria estimada'), ...
    xlabel('Tiempo [seg]'), ylabel('posición x [m]');
subplot(2,2,3), plot(t,x2,t,xh2), legend('velocidad deseada', 'velocidad estimada'), ...
    xlabel('Tiempo [seg]'), ylabel('posición x [m]');
subplot(2,2,4), plot(t,x3,t,xh3), legend('trayectoria deseada', 'trayectoria estimada'), ...
    xlabel('Tiempo [seg]'), ylabel('posición x [m]');
subplot(2,2,2), plot(t,x4,t,xh4), legend('velocidad deseada', 'velocidad estimada'), ...
    xlabel('Tiempo [seg]'), ylabel('posición x [m]');

figure, plot(t,uu);