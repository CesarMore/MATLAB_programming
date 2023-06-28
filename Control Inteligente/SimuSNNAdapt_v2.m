clc
clear all 
close all 

T =0.01;
t = [0:T:30];
x = [0.6 -0.5 1];  %%Condiciones iniciales
%Referencia y derivadas
yr = sin(0.5*t) + 0.5*sin(1.5*t);
yrp = 0.5*cos(t) + 0.75*cos(1.5*t);
yr2p = -0.25*sin(t) - 1.125*sin(1.5*t);
yr3p = -0.0625*cos(t) - 1.6875*cos(1.5*t);
% parametros
% e1 = 0; e2 = 0; e3 = 0;
c = [-1 -1 -1 -1 -1 -1 -1 -1 -1   0    0    0    0    0    0    0    0    0   1  1  1  1  1  1  1  1  1;        
     -1 -1 -1   0    0    0   1  1  1 -1 -1 -1   0    0    0   1  1  1 -1 -1 -1   0    0    0   1  1  1;
     -1   0   1 -1   0   1 -1   0   1 -1   0   1 -1   0   1 -1   0   1 -1   0   1 -1   0   1 -1   0   1];  %vector de centros, numero impar para que sea cero 
% c = [-0.5;0;0.5];
% b = 0.45 ; %bias
% w = rand(9,1); %pesos aleatorios
w = zeros(9,1); %pesos cero
eta = 0.2;   %tasa de entrenamiento 
k1 = 25; k2 = 25; k3 = 25;
K_31 =k1+k2+k3;
K_32 =k1*k2+k1*k3+k2*k3;
K_33 =k1*k2*k3;
K_21 =k1+k2;
K_22 =k2*k1;
Gamma = diag(4);
sigma = 30;%bias
mu = 40;
%Parametros de filtro
psi1 = 1; psi2 = 1; psi3 = 1;
omega1 = 1; omega2 = 1; omega3 = 1;
tau1 = 0.01; tau2 = 0.01; tau3 = 0.01;
T1=psi1*exp(-omega1*t)+tau1;
T2=psi2*exp(-omega2*t)+tau2;
T3=psi3*exp(-omega3*t)+tau3;
%Inicializacion de variable
thetah = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
u(1)=0;

for i=1:length(t)-1
    [tt, xx] = ode45(@modeloRBF3,[t(i) t(i+1)],x(i,:),[],u(i));
    x(i+1,:) = xx(end,:);
    xxx = [x(i+1,1) x(i+1,2) x(i+1,3)]';
    m=27;% No. de nodos
    [a,m]=size(c);% m=27
    for j=1:m
        c_j=c(:,j);
        Xi(j,1)=mu*exp(-(norm(xxx-c_j)^2/(sigma^2)));% Función de activación  mx1
    end    
    z1(i+1) = x(i+1,1) - yr(i);
%     z2(i) = x(i+1,2) - yrp(i);
%     z3(i) = x(i+1,3) - yr2p(i);
    %Ley de control 
    u(i+1) = yr3p(i)-(k1+k2+k3)*(x(i+1,3)+yr2p(i))...
        -(k1*k2+k1*k3+k2*k3)*(x(i+1,2)-yrp(i)) ...
        -(k1*k2*k3)*(x(i+1,1)-yr(i))-(thetah'*Xi);
    %Ley de actualización
    thetah_p = Gamma*((z1(i)*Xi) - (eta*thetah));
    thetah = thetah + (thetah_p*T);
end

figure(1)
plot(t,x(:,1),t,yr)
xlabel('Tiempo (Segundos)'); ylabel('Seguimiento');
legend('x1','yr')
figure(2)
plot(t,z1)
xlabel('Tiempo (Segundos)'); ylabel('error');
legend('z1')
figure(3)
plot(t,u)
xlabel('Tiempo (Segundos)'); ylabel('Control');
legend('u')

function [dxdt] = modeloRBF3(t,x,u)
% --Ejemplo2--
f1=2*x(1)^2*sin(x(2));
f2=x(1)^2+x(1)*x(2)+x(2)*cos(x(1));
f3=x(1)*x(3)+x(2)^2+x(3)*sin(x(2));
dxdt(1,1) = x(2) + f1;
dxdt(2,1) = x(3) + f2;
dxdt(3,1) = u + f3;
end