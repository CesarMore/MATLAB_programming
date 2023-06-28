


A = [-5 4 2 0;
     4 -5 2 0;
     2 2 -8 0];
 R = rref(A)
 
 rrefmovie(A)

clear all 
syms x y 
hold on 
x=-3:0.01:3; 
y1 = 6*x.^2 + 5*x*y - 6*y.^2 + 7; 
plot(x,y1) 
y2=-sqrt(4-4/9*x.^2); 
plot(x,y2) 
grid on
