%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cálculo Simbólico
%
%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formas Canónicas Controlador y Controlable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear
clc
syms s a1 a2 a3 a4 b1 b2 b3 b4
%
Id = eye(4);
%
%------------------------
% Forma canónica controlador
%------------------------
%
% Matrices de la representación de estado
%
Ac = [ 0    1   0   0; 
       0    0   1   0;
       0    0   0   1;
      -a4  -a3 -a2 -a1]
bc = [0; 0; 0; 1]
cc = [b4; b3; b2; b1]
%
%------------------------
% Polos -> valores propios de Ac
%       -> raices del polinomio característico
%------------------------
%
pcc = factor(det(s*Id-Ac)) %polos en forma simbólico
ppcc = charpoly(Ac)        %polos en forma numérico
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MSc = [(s*Id-Ac) bc; -cc.' 0]  %matriz sistema 
pMSc = factor(det(MSc))        %polo de matriz sistema
%
%------------------------
% Matriz de controlabilildad de
% la forma controlador
%------------------------
%
MCc = simplify([bc Ac*bc Ac^2*bc Ac^3*bc])  %matriz de controlabilidad
dMCc = det(MCc)                             %determinante de la matriz de controlabilidad
rank(MCc)                                   %el rango  
%
%------------------------
% Forma canónica controlable
%------------------------
%
% Matrices de la representación de estado
%
Aco = [-a1 1 0 0;
       -a2 0 1 0;
       -a3 0 0 1;
       -a4 0 0 0]
bco = [0; 0; 0; 1]
beta_=inv([1 0 0 0;
                    a1 1 0 0;
                    a2 a1 1 0;
                    a3 a2 a1 1])
beta = simplify(inv([1 0 0 0;
                    a1 1 0 0;
                    a2 a1 1 0;
                    a3 a2 a1 1])*[b1; b2; b3; b4])
cco = [beta(4); beta(3); beta(2); beta(1)]
%
%------------------------
% Polos -> valores propios de Aco
%       -> raices polinomio característico
%------------------------
%
pcco = det(s*Id-Aco)    %polos en forma simbólica
ppcco = charpoly(Aco)   %polos en forma numérica 
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MSco = [(s*Id-Aco) bco; -cco.' 0] % matriz sistema
pMSco = det(MSco)                 % polos de matriz sistema
%
%------------------------
% Matriz de controlabilildad de
% la forma controlable  
%------------------------
%
MCco = [bco Aco*bco Aco^2*bco Aco^3*bco]
dMCco = det(MCco)
rank(MCco)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ejemplo de asignación de polos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear
clc
syms s t
%
Id = eye(3)
vcero = zeros(3,1)
%
%------------------------% Representación de Estado
%------------------------
%
A = [0 1 0;            %como tenemos ceros en la diagonal principal puede ser en la forma
     0 0 1;            %controlador o controlable
     0 0 0]
b = [0; 0; 1]
c = [-1; 0; 1]
x0 = [1; 0; 0]      %condición inicial del sistema 
%
%------------------------
% Funnción de Transferencia Variable de Estado
%------------------------
%
FTve = simplify(c.'*inv(s*Id-A)*b) %
factor(FTve)
[num,den] = ss2tf(A,b,c.',0)
roots(num)
roots(den)
%
%------------------------
% Polos -> valores propios de A
%       -> raices polinomio característico
%------------------------
%
pc = det(s*Id-A) %polinomio forma simbólica
ppc = poly(A)    %polinomio forma numérica
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MS = [(s*Id-A) b; -c.' 0] %matriz sistema
pMS = det(MS)             %polinomio de matriz sistema
factor(pMS)               %factorización del polinomio de matriz sistema
% 
%------------------------
% Matriz de controlabilildad
%------------------------
%
MC = [b A*b A^2*b]      %porque n vale 3 (n=3)
dMC = det(MC)
rank(MC)
%
%------------------------
% Primera proposición:
% FT = ((s-1)(s+1))/((s-1)(s+1)^2) %se hizo esta propuesta (s-1)(s+1)^2
% para eliminar los términos y visualizar que al final en la matriz en la
% forma de Jordan aparece un cero no Hurwitz y eso hace que el sistema no funcione (colapsa) 
%------------------------
%
collect((s+1)^2*(s-1),s) %multiplica el polinomio y simplica, con s variable, simbólica
conv(conv([1 1],[1,1]),[1,-1]) % función convolucion por pares, forma numérica del polonomio
%
% Retroalimentación de estado de
% la primera proposición
%
f1 = [1; 1; -1] % se contruye a partir de la matriz Af, con el 1 de b, en este caso: -(-a3=-1)=1
%                 and so on... s^3 + s^2 - s - 1: {(-a3=-(-1))=1; (-a2=-(-1))=1; (-a3=-(1))=1}
%                 son los coeficientes que se quiere poner 
% Polos lazo cerrado
%
Af1 = A + b*f1.' %matriz en lazo cerrado 
pcLc1 = det(s*Id-Af1) %polinomio característico con matriz en lazo cerrado
factor(pcLc1) 
poly(Af1)
roots(poly(Af1))
%
% Ceros lazo cerrado
%
MSLc1 = [(s*Id-Af1) b; -c.' 0] %matriz sistema en lazo cerrado
pMSLc1 = det(MSLc1)            %polinomio característico del matriz sistema
factor(pMSLc1)
%
% Forma de Jordan del sistema en lazo cerrado
%
[T1, Jf1] = jordan(Af1)
inv(T1)*Af1*T1 - Jf1
bJ1 = inv(T1)*b
cJ1 = T1.'*c
%
% Matriz sistema del sistema en lazo cerrado
% en su forma de Jordan
%
MSLcJ1 = [(s*Id-Jf1) bJ1; -cJ1.' 0]
pretty(MSLcJ1)
%
% Funcióon de transferencia en lazo cerrado
%
FTLc1 = factor(cJ1.'*inv(s*Id-Jf1)*bJ1)
[num1,den1] = ss2tf(Af1,b,c.',0)
roots(num1)
roots(den1)
%
%
%------------------------
% Realizar simulación mdl (Simulink)
%------------------------
%
%
%------------------------
% Segunda proposición:
% FT = ((s-1)(s+1))/(s+1)^3
%------------------------
%
collect((s+1)^3,s)
conv(conv([1 1],[1 1]),[1,1])
%
% Retroalimentación de estado de
% la segunda proposición
%
f2 = [-1; -3; -3]
%AAA
% Polos lazo cerrado
%
Af2 = A + b*f2.'     %matriz en lazo cerrado
pcLc2 = det(s*Id-Af2) %polinomio característico en lazo cerrado
factor(pcLc2)
poly(Af2)
roots(poly(Af2))
%
% Ceros lazo cerrado
%
MSLc2 = [(s*Id-Af2) b; -c.' 0]
pMSLc2 = det(MSLc2)
factor(pMSLc2)
%
% Forma de Jordan del sistema en lazo cerrado
%
[T2, Jf2] = jordan(Af2)
inv(T2)*Af2*T2 - Jf2
bJ2 = inv(T2)*b
cJ2 = T2.'*c
%                   Esto para ver como quedan los modos 
% Matriz sistema del sistema en lazo cerrado en su forma de Jordan
%
MSLcJ2 = [(s*Id-Jf2) bJ2; -cJ2.' 0]
pretty(MSLcJ2)
%
% Funci?n de transferencia en lazo cerrado
%
FTLc2 = cJ2.'*inv(s*Id-Jf2)*bJ2
factor(FTLc2)
[num2,den2] = ss2tf(Af2,b,c.',0)
roots(num2)
roots(den2)
%
%
%------------------------
% Realizar simulaci'on mdl (Simulink)
%------------------------