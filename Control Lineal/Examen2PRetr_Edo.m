%
%
% Examen: parte controlabilidad
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
%------------------
% Representación de Estado
%------------------
%
Ac = [0 1 0;            
      0 0 1;
      0 -1 0]
bc = [0; 0; 1]
cc = [1; 0; -1]
x0 = [-1; 0; 0]     
%
%------------------------
% Función de Transferencia Variable de Estado
%------------------------
%
FTve = simplify(cc.'*inv(s*Id-Ac)*bc) 
factor(FTve)
[num,den] = ss2tf(Ac,bc,cc.',0)
roots(num)
roots(den)
%
%------------------------
% Polos -> valores propios de A
%       -> raices polinomio característico
%------------------------
%
pc = det(s*Id-Ac)
ppc = charpoly(Ac)
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MS = [(s*Id-Ac) bc; -cc.' 0] 
pMS = det(MS)             
factor(pMS)               
% 
%------------------------
% Matriz de controlabilildad
%------------------------
%
MCc = [bc Ac*bc Ac^2*bc]      
dMCc = det(MCc)
rank(MCc)
%-----------------------------------------------------------------------------
% Forma canonica de controlabilildad
%---------------------------------------------------------------------------------
Aco = [0 1 0;            
      -1 0 1;
       0 0 0]
bco = [ 0; 0; 1]
beta = (inv([1 0 0;
             0 1 0 ;
             1 0 1])*[-1; 0; 1])
cco = [beta(3); beta(2); beta(1)]

FTve = simplify(cco.'*inv(s*Id-Aco)*bco) 
factor(FTve)
[num,den] = ss2tf(Aco,bco,cco.',0)
roots(num)
roots(den)
%
%------------------------
% Polos -> valores propios de A
%       -> raices polinomio característico
%------------------------
%
pc = det(s*Id-Aco) 
ppc = charpoly(Aco)   
%
%------------------------
% Ceros -> raices del determinante de
%          la Matriz Sistema
%------------------------
%
MS = [(s*Id-Aco) bco; -cco.' 0] 
pMS = det(MS)             
factor(pMS)

%------------------------
% Matriz de controlabilildad
%------------------------
%
MCco = [bco Aco*bco Aco^2*bco]      
dMCco = det(MCco)
rank(MCco)

collect((1+s)*(s+1)^2,s)
%conv(conv([1 1],[1 1]),[1,1])
% Retroalimentación de estado de
% la primera proposición
%
f1 = [-1; -2;-3] 
%
% Polos lazo cerrado
%
Af1 = Ac + bc*f1.'  
pcLc1 = det(s*Id-Af1) 
factor(pcLc1) 
poly(Af1)
roots(poly(Af1))
%
% Ceros lazo cerrado
%
MSLc = [(s*Id-Af1) bc; -cc.' 0] 
pMSLc = det(MSLc)            
factor(pMSLc)
%
% Forma de Jordan del sistema en lazo cerrado
%
[T1, Jf1] = jordan(Af1)
inv(T1)*Af1*T1 - Jf1
bJ1 = inv(T1)*bc
cJ1 = T1.'*cc
%
% Matriz sistema del sistema en lazo cerrado
% en su forma de Jordan
%
MSLcJ1 = [(s*Id-Jf1) bJ1; -cJ1.' 0]
pretty(MSLcJ1)
%
% Función de transferencia en lazo cerrado
%
FTLc1 = factor(cJ1.'*inv(s*Id-Jf1)*bJ1)
[num1,den1] = ss2tf(Af1,bc,cc.',0)
roots(num1)
roots(den1)
%
%
%------------------------
% Realizar simulación mdl (Simulink)
%------------------------
%
%
%
% Examen: parte observabilidad
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear
clc
syms s
%
Id = eye(3);
vcero = zeros(3,1);
%
%------------------------
% Representaci'on de Estado
%------------------------
%
Ac = [0 1 0;          
     0 0 1;
     0 -1 0]
bc = [0; 0; 1]
cc = [1; 0; -1]
x0 = [-1; 0; 0]
%
%------------------------
% Polos -> valores propios de A
%       -> raices polinomio característico
%------------------------
%
pc = det(s*Id-Ac)
ppc = charpoly(Ac)
%
%------------------------
% Ceros -> raices del determinante de la Matriz Sistema
%------------------------
%
MS = [(s*Id-Ac) bc; -cc.' 0]
pMS = det(MS)
%
%
%------------------------
% Matriz de controlabilidad cotrolador 
%------------------------
MCc = [bc Ac*bc Ac^2*bc]      
pMCc = det(MCc)
dMCc = charpoly(MCc)
rank(MCc)
%------------------------
% Matriz de observador controlador 
%------------------------
%
MOc = [cc.'; cc.'*Ac; cc.'*Ac^2]
dMOc = det(MOc)
pMoc= charpoly(MOc)
rank(MOc)
%
%
%------------------------
% Cálculo de la retroalimentación de estado
%------------------------
collect((1+s)*(s+1)^2,s)
%

%------------------------
%
% Transformación a la forma observador
%
Ao = [0 0 0;
      1 0 -1;
      0 1 0]
bo = [1; 0; -1]
co = [0; 0; 1]
[Ao; co.']
%
factor(co.'*inv(s*Id-Ao)*bo)
%
%
MOo = [co.'; co.'*Ao; co.'*Ao^2]
dMOo = det(MOo)
pMCc = charpoly(MOo)
rank(MOo)
%
%
% Transformación a la forma controlable
%
%
%
Aco = [0 1 0;            
      -1 0 1;
       0 0 0]
bco = [ 0; 0; 1]
beta = (inv([1 0 0;
             0 1 0 ;
             1 0 1])*[-1; 0; 1])
cco = [beta(3); beta(2); beta(1)]
%
%
% Transformación a la forma observable
%
%
Aob = [0 1 0;            
      -1 0 1;
       0 0 0]
bob = [ 0; 0; 1]
beta = (inv([1 0 0;
             0 1 0 ;
             1 0 1])*[-1; 0; 1])
cob = [beta(3); beta(2); beta(1)]
%
%
%Funcion de tranferencia--------------------------------------------------
%
% 
FTve = simplify(cco.'*inv(s*Id-Aco)*bco) 
factor(FTve)
[num,den] = ss2tf(Aco,bco,cco.',0)
roots(num)
roots(den)
%
%
MCco = [bco Aco*bco Aco^2*bco]      
%
MOco = [cco.';cco.'*Aco;cco.'*Aco^2]
dMOco = det(MOco)
rank(MOco)
%
%
% Retroalimentación de estado--------------------------------------------
%
f = [-1; -2; -3]
%
Afc = Ac + bc*f.'
FTve = simplify(cc.'*inv(s*Id-Ac)*bc) 
factor(FTve)
factor(det(s*Id-Afc))
poly(Afc)
eig(Afc)
%
%
% Cambio de base de controlador a controlabilidad--------------------------
%
%
T1 = MCc*inv(MCco)
inv(T1)*Ac*T1 - Aco  
inv(T1)*bc - bco     
cc.'*T1 - cco.'      
%
%%retro en controlabilidad-------------------------------------------------
%
fco=T1*f    
%
Afco = Aco + bco*fco.'
factor(FTve)
factor(det(s*Id-Afco))
poly(Afco)
eig(Afco)
%
%------------------------
% Inyección de salida
%------------------------
collect((5+s)*(s+5)^2,s)
%
ko = [-125; -74; -15]
%
Ako = Ao + ko*co.'
factor(det(s*Id-Ako))
%-----------------------------------------------
%CAMBIO DE BASE DE OBSERVABOR A CONTROLABLE
T2 = inv(MOco)*MOo 
%
%inyeccion de salida controlable
kco = T2*ko
%
%
%inyeccion de salida controlador
T3 = inv(MOc)*MOo  % base de 
kc = T3*ko
%
 
%
%
%
Akc = Ac + kc*cc.'
factor(det(s*Id-Akc))
%
%
Akco = Aco + kco*cco.'          
factor(det(s*Id-Akco))
%
%------------------------
% Sistema en Lazo cerrado
%------------------------
Alcc = [ Ac   bc*f.' ;
       -kc*cc.' (Afc+kc*cc.')]
blcc = [bc; bc]
clcc = [cc; vcero]
%
%
Alcco = [ Aco   bco*fco.' ;
         -kco*cco.' (Afco+kco*cco.')]
blcco = [bco; bco]
clcco = [cco; vcero]
%------------------------
% Polos -> valores propios de Alc
%       -> raices polinomio característico
%------------------------
%
I6 = eye(6);
pclcc = det(s*I6-Alcc)
factor(pclcc)
%
%
pclcco = det(s*I6-Alcco)
factor(pclcco)
%
%------------------------------------------------------
% Ceros -> raices del determinante de la Matriz Sistema
%------------------------------------------------------
%
MSlcc = [(s*I6-Alcc) blcc; -clcc.' 0]
pMSlcc = det(MSlcc)
factor(pMSlcc)
%
%
MSlcco = [(s*I6-Alcco) blcco; -clcco.' 0]
pMSlcco = det(MSlcco)
factor(pMSlcco)
%
%---- --------------------
% Matriz de controlabilildad
%------------------------
%
MClc = [blcc Alcc*blcc Alcc^2*blcc Alcc^3*blcc Alcc^4*blcc Alcc^5*blcc]
dMClc = det(MClc)
rkMClc = rank(MClc)
factor(pclcc)
factor(pMSlcc)
%
%
MClco = [blcco Alcco*blcco Alcco^2*blcco Alcco^3*blcco Alcco^4*blcco Alcco^5*blcco]
dMClco = det(MClco)
rkMClco = rank(MClco)
factor(pclcco)
factor(pMSlcco)
%
%------------------------
% Matriz de observabilildad
%------------------------
%
MOlcc = [clcc.';
        clcc.'*Alcc;
        clcc.'*Alcc^2;
        clcc.'*Alcc^3;
        clcc.'*Alcc^4;
        clcc.'*Alcc^5]
dMOlcc = det(MOlcc)
rkMOlcc = rank(MOlcc)
%factor(pMOlcc)
%
%
MOlcco = [clcco.';
        clcco.'*Alcco;
        clcco.'*Alcco^2;
        clcco.'*Alcco^3;
        clcco.'*Alcco^4;
        clcco.'*Alcco^5]
dMOlcco = det(MOlcco)
rkMOlcco = rank(MOlcco)
factor(pclcco)
factor(pMSlcco)
%
%--------------------------------------------
% Forma de Jordan del Sistema en Lazo cerrado
%--------------------------------------------
%
[T4, Jlcc] = jordan(Alcc)
round(inv(T4)*Alcc*T4 - Jlcc)
bJlcc = inv(T4)*blcc
cJlcc = T4.'*clcc
MSlcc =[(s*I6-Jlcc) bJlcc;
        -cJlcc.'    0  ]
vpa(MSlcc,2)
round(MSlcc)
%
%
[T5, Jlcco] = jordan(Alcco)
round(inv(T5)*Alcco*T5 - Jlcco)
bJlcco = inv(T5)*blcco
cJlcco = T5.'*clcco
MSlcco =[(s*I6-Jlcco) bJlcco;
        -cJlcco.'    0  ]
vpa(MSlcco,2)
round(MSlcco)
%
% 
