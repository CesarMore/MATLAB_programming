%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cálculo Simbólico
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Forma canónica observador
%------------------------
%
% Matrices de la representación de estado
%
Ao = [0 0 0 -a4;
      1 0 0 -a3;
      0 1 0 -a2;
      0 0 1 -a1]
bo = [b4; b3; b2; b1]
co = [0; 0; 0; 1]
%
%------------------------
% Polos -> valores propios de Ao
%       -> raices del polinomio característico
%------------------------
%
pco = factor(det(s*Id-Ao)) %polos en forma simbólico
ppco = charpoly(Ao)        %polos en forma numérico
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MSo = [(s*Id-Ao) bo; -co.' 0]  %matriz sistema 
pMSo = factor(det(MSo))        %polo de matriz sistema
%
%------------------------
% Matriz de observabilidad de
% la forma observador
%------------------------
%
MCo = simplify([co.'; co.'*Ao; co.'*Ao^2; co.'*Ao^3])  %matriz de controlabilidad
dMCo = det(MCo)                             %determinante de la matriz de controlabilidad
rank(MCo)                                   %el rango  
%
%------------------------
% Forma canónica observable
%------------------------
%
% Matrices de la representación de estado
%
Aob = [-a1 -a2 -a3 -a4;
        1   0   0   0;
        0   1   0   0;
        0   0   1   0]
cob = [0; 0; 0; 1]
beta = simplify(inv([1 0 0 0;
                    a1 1 0 0;
                    a2 a1 1 0;
                    a3 a2 a1 1])*[b1; b2; b3; b4])
bob = [beta(4); beta(3); beta(2); beta(1)]
%
%------------------------
% Polos -> valores propios de Aob
%       -> raices polinomio característico
%------------------------
%
pcob = det(s*Id-Aob)    %polos en forma simbólica   
ppcob = charpoly(Aob)   %polos en forma numérica
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MSob = [(s*Id-Aob) bob; -cob.' 0] % matriz sistema
pMSob = det(MSob)                 % polos de matriz sistema
%
%------------------------
% Matriz de observabilidad de
% la forma observable  
%------------------------
%
MOob = simplify([cob.'; cob.'*Aob; cob.'*Aob^2; cob.'*Aob^3])
dMOob = det(MOob)
rank(MOob)
%