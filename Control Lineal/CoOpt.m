%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ejemplo Control Óptimo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all
clc
syms s
Id = eye(3)
vcero = zeros(3,1)
%
%------------------------
% Representación de Estado
%------------------------
%
A = [0 1 0;
    0 0 1;
    0 0 0]
b = [0; 0; 1]
c = [-1; 0; 1]
x0 = [1; 0; 0]
%
%------------------------
% Polos -> valores propios de A
%------------------------
%
pc = det(s*Id-A)
vpA = eig(A)
%
%------------------------
% Ceros -> raices del determinante de la Matriz Sistema
%------------------------
%
MS = [(s*Id-A) b; -c.' 0]
pMS = (det(MS))
factor(pMS)
%
%------------------------
% Función de Transferencia Variable de Estado
%------------------------
%
FTve = simplify(c.'*inv(s*Id-A)*b)
%
%------------------------
% Matriz de observabilildad
%------------------------
%
MO = [c.'; c.'*A; c.'*A^2]
dMO = det(MO)
%
%------------------------
% Matriz de controlabilildad
%------------------------
%
MC = [b A*b A^2*b]
dMC = det(MC)
%
%%
%------------------------
% Proposición de asignación de polos:
% FT = ((s-1)(s+1))/(s+1)^3
% FT = ((s-1)(s+1))/(s+10)^3
% FT = ((s-1)(s+1))/(s+0.1)^3
%------------------------
%
collect((s+1)^3,s)
%collect((s+10)^3,s)
%collect((s+0.1)^3,s)
%
% Retroalimentación de estado asignación de polos
%
%fp = [-1; -3; -3]
%fp = [-1000; -300; -30]
fp = [-1/1000; -3/100; -3/10]
%
% Polos lazo cerrado
%
Afp = (A + b*fp.')
vpAfp = eig(Afp)     %vp valores propios
%
% Ceros lazo cerrado
%
MSLcp = [(s*Id-Afp) b; -c.' 0]
pMSLcp = (det(MSLcp))
factor(pMSLcp)
%
% Función de transferencia en lazo cerrado
% 
SisLRp= ss(Afp,b,c.',0)
TFLRrp = zpk(SisLRp)
%
%Esto se obtuvo con la retro de estado, 
%posicionando los polos en s+1
%------------------------
% 
%%
% Ahora comparación con Control Óptimo
%------------------------
% Proposición <Control Óptimo>:
% Índice de desempeño o el criterio Dinámica salida-entrada
% J = (1/2)\int_{0}^{\infty}(y^2 + \rho u^2)dt
%------------------------
%
Q1e2 = (c.')       %Q1e2 es la Q 1/2
Q = (Q1e2.'*Q1e2)  % con esto se minimiza el criterio
                   % que es lo mismo c.'*c     
vpQ = eig(Q)
%r = 1 
%r = 10
r = 1/10
%
%------------------------
% El par (A,b) es controlable
rank([b,A*b,A^2*b])  %[b,A*b,A^2*b] matriz de controlabilidad
McPBH = [(s*Id-A),b] %matriz de controlabilidad de prueba PBH
Tc = [1,   0,  0, 0;
      s,  -1,  0, 0;
     s^2, -s, -1, 0;
    -s^3, s^2, s, 1]
det(Tc)
McPBH2 = McPBH*Tc 
%------------------------
% El par (Q1/2,A) es observable
rank([Q1e2;Q1e2*A;Q1e2*A^2])
MoPBH = [Q1e2;(s*Id-A)]
To = [-1,  0,   0,  0;
     -s,  -1,   0,  1;
     -s^2,-s,  -1,  s;
     -s^3,-s^2,-s,-1+s^2]
det(To)
MoPBH2 = simplify(To*MoPBH)
%
%------------------------
% Ecuación de Riccati:
% A.'*P + P*A - P*b*r^(-1)*b.'*P + Q = 0
%
% Retroalimentación Óptima: fo.' = r^(-1)*b.'*P
%------------------------
%
[P,L,G] = care(A,b,Q,r,vcero,Id)  %G es lavretroalimentación
P                                 %P es la solucion de la matriz P
vpP = eig(P)                      %L es el espectro (valores propios)  
% P simetrica positiva definida;   de la matriz retroalimentado por G
%------------------------
% Retroalimentación de estado control Óptimo%
G                  %1)las 3 retros son los mismo
r^(-1)*b.'*P       %2)
K = lqr(A,b,Q,r)   %3)
fo = G.' 
%
%------------------------
% Polos lazo cerrado
%------------------------
%
Afo = A - b*fo.'
pcLco = vpa(det(s*Id-Afo),4)
poly(Afo)
vpAfo = eig(Afo)
%
%------------------------
% Ceros lazo cerrado
%------------------------
%
MSLco = vpa([(s*Id-Afo) b; -c.' 0],4)
pMSLco = det(MSLco)
factor(pMSLco)
%
%
%------------------------
% Función de transferencia lazo cerrado
%------------------------
%
SisLRo= ss(Afo,b,c.',0)
TFLRro = zpk(SisLRo)
[knumLC,denLC] = ss2tf(Afo,b,c.',0)
ceros = roots(knumLC)
polos = roots(denLC)
%
%------------------------
% Comparación polos en lazo cerrado
%------------------------
%
vpAfp
r
vpAfo
[vpAfp,vpAfo]
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realizar simulación mdl (Simulink)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
