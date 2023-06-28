clear all
clc
close all

%
C = 10;
a = 1.1*C;
l = sqrt(C);

T = 0.01;
t = [0:T:20];

f = 3*cos(t);
df = -3*sin(t);

x = [0];
u1 = [0];
u = [0];
u1p = [0];

%Integrando para obtener u1 usando métdo de Euler
for i=1:length(t)-1

    u1p(i+1) = -a*sign(x(i)-f(i));
    u1(i+1) = u1(i) + T*(u1p(i));
    
    u(i+1) = u1(i+1) - l*sqrt(abs(x(i)-f(i)))*sign(x(i)-f(i));
    x(i+1) = x(i) + T*(u(i+1));
    
end

figure(1)
plot(t,u,t,df)
legend('Diferenciador Levant','Derivada exacta de f')

















