clc
clear all 
close all 



%[tt, xx] = ode45(@modeloPVTOL, [0:0.1:5], [1 0.1]); 
[tt, xx] = ode45(@modeloPVTOL, [0:0.1:5], [1 0.1]); 
%[tt, xx] = ode45(@modeloPVTOL, [0:0.1:5], [1 0], [], 0.1); 


plot(tt,xx(:,1),tt,xx(:,2))
grid
legend('x_1','x_2')
%%
%Gráfica de la solución del vector de estado 
%
%dx(t)/dt =A*x(t) + b*u(t)
%
%x(t) = e^(A*t).*x(0) + integral(e^(A(t-tau))*b*u(tau)dtau
%
clc
clear all
close all

t = 0:0.1:5
x1= t
plot(t,x1)
hold on
x2 = 0.1*t
plot(t, x2)
legend('x_1','x_2')









%%
clear all 
close all 
clc

syms xx1 xx2 

xi = [1;0.1]

for t=0:0.1:5
    exp = [1 t;
           0 1]
    x= exp*xi
    xx1 = [xx1 ; x(1)]
    xx2 = [xx2 ; x(2)]
    %tt = [tt ; t]
end


plot(t,xx1,t,xx2)
legend('x_1','x_2')
    


