%%
%Gráfica de la solución del vector de estado 
%
%dx(t)/dt =A*x(t) + b*u(t)
%
%x(t) = e^(A*t).*x(0) + integral(e^(A(t-tau))*b*u(tau)dtau
%

close all
clear all 
clc
syms s z

tt = [0]

vecx1 = {}
vecx2 = {}
for T = 0:0.5:5
    
    u = 0.1
    A = [1 0;0 1]
    x0 = [1;0.1]
    B = [0;1]

    inversa = s*eye(2)-A
    X = inv(inversa)*B*u
    XL = ilaplace(X,z)
    XLs = subs(XL,z,T)

    L = ilaplace(inv(inversa),z)
    Ls = subs(L,z,T)

    est = Ls*x0+XLs
    
    vecx1 = [vecx1;est(1)]
    vecx2 = [vecx2;est(2)]
    
    tt = [tt;T]
end

x1e = vecx1
x2e = vecx2

tt = [tt(2);tt(3);tt(4);tt(5);tt(6);tt(7);tt(8);tt(9);tt(10);tt(11);tt(12)]
%tt = [0] 
hold on
plot(tt,x1e,tt,x2e)
grid
