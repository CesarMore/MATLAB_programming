clc
clear all 
close all 

[tt, xx] = ode45(@modeloPVTOL_1, [0:0.1:5], [0 0 0 0 0 0],[], 0 ,0 ); 


plot(tt,xx(:,1),tt,xx(:,2),tt,xx(3),tt,xx(4),tt,xx(5),tt,xx(6))
grid
legend('x_1','x_2','x_3','x_4','x_5','x_5')