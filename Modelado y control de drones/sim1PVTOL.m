clc 
clear all 
%close all

T = 0.1;              % periodo de simulación 
t = [0:T:5];          % tiempo de simulación
x = [0 0 0 0 0 0];
zd = 2;
k1z = 2;
k2z = 1;
g = 9.81;
m = 1;

% 1.- Medir las variables a controlar 
[tt, xx] = ode45(@modelo1PVTOL,[0 T],x(1,:),[],0,0);
x(2,:) = xx(end,:)

% 2.- Calcular la ley de control 
ez = x(2,2) - zd;
ezp = x(2,5);
fq = (m/cos(x(2,1)))*( -k1z*ezp - k2z*ez +g )

% 3.- Enviar la ley de contro a los actuadores

[tt, xx] = ode45(@modelo1PVTOL,[0 T],x(2,:),[],fq,0);
x(3,:) = xx(end,:);










