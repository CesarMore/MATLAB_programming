clc
clear all
close all

%t = 30;
t = 1*[-pi/2:0.01:pi/2];
y = 2*cos(t).*sin(2*t);    %función a aproximar
c = (pi/2)*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];  %vector de centros, numero impar para que sea cero 
b = 0.45 ; %bias
w = rand(9,1);
eta = 0.1;   %tasa de entrenamiento 
inc = 1;

for k=1:30  % k numero de etapas
    for i=1:length(t)
        for j=1:length(c)
            h(j,1) = exp(- norm(t(i)- c(j))^2/ (2*b^2));
%            f1(i) = exp(-norm(x(i) -c(1))^2 / (2*b^2));
%            f2(i) = exp(-norm(x(i) -c(1))^2 / (2*b^2));
%            f3(i) = exp(-norm(x(i) -c(1))^2 / (2*b^2));
%            f4(i) = exp(-norm(x(i) -c(1))^2 / (2*b^2));
%            f5(i) = exp(-norm(x(i) -c(1))^2 / (2*b^2));
        end
    yNN(i) = w'*h;

    w = w + eta*(y(i) - yNN(i))*h;
    W(inc,:) = w;
    inc = inc+1;
    end
    E(k) = inv(2)*norm(y-yNN)^2; %error de etapa
end

figure(1)
plot(t, y, t, yNN)
legend('y', 'yNN')
figure(2)
plot(W)
title('Pesos')
figure(3)
plot(E)
title('Error de entrenamiento')
