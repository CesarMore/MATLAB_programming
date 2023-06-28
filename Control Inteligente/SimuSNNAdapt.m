clc
clear all 
close all 

T =0.01;
t = [0:T:10];
x = [0.5 -0.3 0.4];  %%Condiciones iniciales

%Referencia y derivadas
yr = 0.5*sin(t);
yrp = 0.5*cos(t);
yr2p = -0.5*sin(t);
yr3p = -0.5*cos(t);

% parametros
% e1 = 0; e2 = 0; e3 = 0;
c = [-0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5   0    0    0    0    0    0    0    0    0   0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5;        
     -0.5 -0.5 -0.5   0    0    0   0.5  0.5  0.5 -0.5 -0.5 -0.5   0    0    0   0.5  0.5  0.5 -0.5 -0.5 -0.5   0    0    0   0.5  0.5  0.5;
     -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5 -0.5   0   0.5];  %vector de centros, numero impar para que sea cero 
%c = [-0.5;0;0.5];
% b = 0.45 ; %bias
% w = rand(9,1); %pesos aleatorios
w = zeros(9,1); %pesos cero
eta = 0.2;   %tasa de entrenamiento 
k1=12; k2=10; k3=12;
K_31=k1+k2+k3;
K_32=k1*k2+k1*k3+k2*k3;
K_33=k1*k2*k3;
K_21=k1+k2;
K_22=k2*k1;
Gamma = diag(4);
sigma = 10;
mu = 10;

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

%     T1=psi1*exp(-omega1*t)+tau1;
%     T2=psi2*exp(-omega2*t)+tau2;
%     T3=psi3*exp(-omega3*t)+tau3;
%     Yrsp=zeros(3,1);
%     Yrsp(1,1)=(yr-Yr)/T1;
%     Yrs(i+1) =  Yrs(i) + Yrsp(i)+T;
%     Yrsp(2,1)=(yrp-Yrp)/T2;
%     
%     Yrsp(3,1)=(yr2p-Yr2p)/T3;
%     Yrsp(4,1)=(yr3p-Yr3p)/T4;
%     Yr=Yrs(1,1);
%     Yrp=Yrs(2,1); 
%     Yr2p=Yrs(3,1);
%     %
%     z1=x1-Yr;%1x1
%     xtilde=x1-Yr;
%     f1=x1^3;
%     x1p=x2+f1;
%     xtildep=x1p-yrp;
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
% --Ejemplo1--
f1 = x(1)^3;
f2 = x(1)^2 + x(2)^2;
f3 = 0;
dxdt(1,1) = x(2) + f1;
dxdt(2,1) = x(3) + f2;
dxdt(3,1) = u + f3;
end