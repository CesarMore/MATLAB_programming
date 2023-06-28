clc
clear all
close all

t = 1*[-pi/2:0.01:pi/2];
y = 2*cos(t).*sin(2*t);    %función a aproximar
c = (pi/2)*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];  %vector de centros, numero impar para que sea cero 
b = 0.45 ; %bias
w = rand(9,1);
eta = 0.1;   %tasa de entrenamiento 
inc = 1;

t1 = [-2:0.01:2];
y1 = sin(2*t1); 

for i=1:length(t1)
    for j=1:length(c)
        h(j,1) = exp(- norm(t1(i)- c(j))^2/ (2*b^2));
    end
    yNN1(i) = w'*h;
end

figure(1)
plot(t1, y1, t1, yNN1)
legend('y1', 'yNN1')

figure(2)
plot(w)
title('Pesos')
figure(3)
plot(c)
title('Error de entrenamiento')