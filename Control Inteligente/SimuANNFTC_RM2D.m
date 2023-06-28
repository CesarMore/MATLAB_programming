clc
clear all 
close all 

T =0.01;
t = [0:T:30];
x = [0.6 -0.5 1];  %%Condiciones iniciales
%Referencia y derivadas
q1d = sin(2*pi*t);
q2d = sin(2*pi*t);
q1dp = 2*pi*cos(2*pi*t);
q2dp = 2*pi*cos(2*pi*t);
% parametros
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
alpha1 = 0.7;
alpha2 = 0.8235;
Kp = [100 0; 0 100];
Kd = [500 0; 0 500];
fv1 = 0.534684;
fv2 = 0.001;
fc1 = 0.81958;
fc2 = 0.002;
%Inicializacion de variable
thetah = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
u(1)=0;

for i=1:length(t)-1
    [tt, xx] = ode45(@modeloANNFTC,[t(i) t(i+1)],x(i,:),[],u(i));
    x(i+1,:) = xx(end,:);
    xxx = [x(i+1,1) x(i+1,2) x(i+1,3)]';
    m=27;% No. de nodos
    [a,m]=size(c);% m=27
    for j=1:m
        c_j=c(:,j);
        Phi(j,1) = exp(-(norm(xxx - c_j)^2/(sigma^2)));% Función de activación  mx1
    end    
    z1(i+1) = x(i+1,1) - yr(i);

    %Ley de control 
    error1 = x(i+1,1)-q1d;
    error2 = x(i+1,3)-q2d;
    errorp1 = x(i+1,2)-q1dp; 
    errorp2 = x(i+1,4)-q2dp;
    Sig_e = [(abs(error1))^alpha1*sgn(error1); (abs(error2))^alpha1*sgn(error2)];
    Sig_ep = [(abs(errorp1))^alpha2*sgn(errorp1); (abs(errorp2))^alpha2*sgn(errorp2)];
    p1 = 0.0398; p2 = 0.0026; p3 = -0.0015; p4 = 0.0081;
    m11 = p1 + p2 + 2*p3*cos(x(i+1,3)) - 2*p4*sin(x(i+1,3));
    m12 = p2 + p3*cos(x(i+1,3)) - p4*sin(x(i+1,3));
    m22 = p2;
    M = [m11 m12; m12 m22];
    C = [-2*b*x(i+1,4) -b*x(i+1,4); b*x(i+1,2) 0];
    qp = [x(i+1,2); x(i+1,4)];
    qd2p = [q1dp;q2dp];
    b = p3*sin(x(i+1,3)) + p4*cos(x(i+1,3));
    f1 = fv1*x(i+1,2) + fc1*sgn(x(i+1,2));
    f2 = fv2*x(i+1,4) + fc2*sgn(x(i+1,4));
    F = [f1 ; f2];
    tau = M*( qd2p - Kp*(Sig_e)^alpha1 - Kd*(Sig_ep)^alpha2) + C*qp + F;
    tau1 = tau(1,1);
    tau2 = tau(2,1);
    %Ley de actualización
    thetah_p = Gamma*(Phi*z1);
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

