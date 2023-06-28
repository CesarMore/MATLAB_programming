%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Programa MatLab correspondiente a
%
% a Asignación de Polos:
%
% Teoría de Anillos
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
clc
syms s a e
N = 2*(s+a)
M = collect((s+1)*(s^2+s+1),s)
Q = collect((s+1)*(s+1/2)^2,s) % Q polinomio Hurwitz dado
%
%%%%%%%%%%%%%%%%%%%%%
%%
%
%Algoritmo de división de Euclides
%
%S = 1
%R = Q-M*S
[S,R] = quorem(Q,M,s)       %Q = S*M+R 
factor(R)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
% Ecuacion Diofantina
%
a1=2; a2=2; a3=1;
b0=0; b1=0; b2=2; b3=2*a;
X = [1 0 0 b0 0 0;
    a1 1 0 b1 b0 0;
    a2 a1 1 b2 b1 b0;
    a3 a2 a1 b3 b2 b1;
    0 a3 a2 0 b3 b2;
    0 0 a3 0 0 b3]
factor(det(X))

Xinv = simplify(inv(X))
pretty(Xinv)
pLC = collect((s+1/e)^2*Q,s)  % polinomio, 
q = [1;                       % vector q de la matriz Silvester  
(2/e + 2);
(4/e + 1/e^2 + 5/4);
(5/(2*e) + 2/e^2 + 1/4);
(1/(2*e) + 5/(4*e^2));
1/(4*e^2)]
coef = simplify(inv(X)*q)     % se puede considerar el vector v de la matriz Silvester
Sd = collect(coef(1)*s^2 + coef(2)*s + coef(3),s) %S Diofantina para construir matriz Silvester
Rd = collect(coef(4)*s^2 + coef(5)*s + coef(6),s) %R Diofantina para construir matriz Silvester
Qob = simplify(Sd*M+Rd*N)                         % Q observador o Diofantina (Q=S*M+R*N)
factor(Qob) % Importante verificar que factor(Qob)
factor(pLC) %  y factor(pLC) sean iguales. 
FT = simplify((1/Sd)*((N/M)/(1+(Rd/Sd)*(N/M))))
pretty(FT)
% FT
% Polinomios para a = 1/2
%
Sfm = collect(subs(Sd,[a e],[1/2 0.1]),s)
Rfm = collect(subs(Rd,[a e],[1/2 0.1]),s)
Nfm = subs(N,a,1/2)
FTfm = simplify((1/Sfm)*((Nfm/M)/(1+(Rfm/Sfm)*(Nfm/M))))
factor(FTfm)
% FT
% Polinomios para a = -1/2
%
Sfnm = collect(subs(Sd,[a e],[-1/2 0.1]),s)
Rfnm = collect(subs(Rd,[a e],[-1/2 0.1]),s)
Nfnm = subs(N,a,-1/2)
FTfnm = simplify((1/Sfnm)*((Nfnm/M)/(1+(Rfnm/Sfnm)*(Nfnm/M))))
factor(FTfnm)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Resultados del Programa
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Datos de Simulación MatLab
%
% correspondiente a
%
% Asignación de Polos:
%
% Teoría de Anillos
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
clc
%a = 1/2
a = -1/2
%
e = 1/10
N = [2 2*a]
M = conv([1 1],[1 1 1])
Qe = conv([1 1],conv([1 1/2],[1 1/2]))
%
% Algoritmo de divisi'on de Euclidez
%
%Se = 1
%Re = Qe-conv(M,Se)

[Se,Re] = deconv(Qe,M)
De = conv(Se,N) 
%
% Ecuacion Diofantina
%
a1=2; a2=2; a3=1;
b0=0; b1=0; b2=2; b3=2*a;
X = [1 0 0 b0 0 0;
    a1 1 0 b1 b0 0;
    a2 a1 1 b2 b1 b0;
    a3 a2 a1 b3 b2 b1;
    0 a3 a2 0 b3 b2;
    0 0 a3 0 0 b3]
det(X)
Xinv = inv(X)
Qd = conv(conv([1 1/e],[1 1/e]),Qe)
q = [1;
(2/e + 2);
(4/e + 1/e^2 + 5/4);
(5/(2*e) + 2/e^2 + 1/4);
(1/(2*e) + 5/(4*e^2));
1/(4*e^2)]
coef = Xinv*q
Sd = [coef(1) coef(2) coef(3)]
Rd = [coef(4) coef(5) coef(6)]
num = N
roots(N)
den = conv(M,Sd) + conv([0 0 1],conv(N,Rd))
roots(den)
gg = num(2)/den(6)

