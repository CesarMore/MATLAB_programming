% Graficar 3 funciones a la vez para observar sus acotamiento
%fplot(@(x) funcion)
%Funciones : x,2x,3x
fplot(@(x) sin(x))
hold on
grid minor
fplot(@(x) 1)
legend
hold on
fplot(@(x) -1)
legend
hold on
fplot(@(x) 0,'k') %Marcar el eje X
hold on
plot([0,0],[10,-10],'k')%Marcar el eje Y
hold off

%Matriz 2x2
A=[ 5 0 ;
    0  1 ]
Valores_propios=eig(A)


