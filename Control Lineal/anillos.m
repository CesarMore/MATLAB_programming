clear
clc
syms s a e


e = 1/5
N = [0 -4 0 1]
M = [1 0 -1 0]
Qcoq = collect((s+1)^3,s)
Qe = [1 3 3 1]

roots(N)
roots(M)

% 
% Algoritmo de divisi'on de Euclidez
%
%Se = 1
%Re = Qe-conv(M,Se)
[Se,Re] = deconv(Qe,M)
comp = Se*M + Re
De = conv(Se,N)
entrada = conv(Qe,N)
%%%%Comprobando la TF lc
%N = collect(N(2)*s^2 + N(4),s)
%M = collect(M(1)*s^3 + M(3)*s,s)
%Se = 1
%Re = collect(Re(2)*s^2 + Re(3)*s+Re(4),s)
%FT = simplify((1/Se*N)((N/M)/(1+(Re/Se*N)(N/M))))


% Ecuacion Diofantina
%
a1=0; a2=-1; a3=0;
b0=0; b1=-4; b2=0; b3=1;
X = [1 0 0 b0 0 0;
a1 1 0 b1 b0 0;
a2 a1 1 b2 b1 b0;
a3 a2 a1 b3 b2 b1;
0 a3 a2 0 b3 b2;
0 0 a3 0 0 b3]
det(X)
Xinv = inv(X)
Y = round(Xinv,2)

%Qd = conv(conv([1 1/e],[1 1/e]),Qe)
%Qd = conv((conv(conv([1 1/e],[1 1/e]),[2 1])/(conv([0 0 1],[1 1]))),Qe)

Qd = conv(conv([1 2.5],[1 2.5]),conv(conv([1 1],[1 1]),[2 1]))
%Qd2 = conv(conv([1 5],(conv([1 5],[1 5]))),conv([1 1],[1 1]))
%Qd = conv(conv(.5,[2 1]),conv(conv([1 5],[1 5]),conv([1 1],[1 1]))) 

%q = [2  ;  25  ; 104  ; 166  ; 110  ;  25]
q = [ 2.0000  ; 15.0000 ; 41.5000 ;  52.2500  ; 30.0000  ;  6.2500]
%q = q/2

coef = Xinv*q
%%coef = Y*q2

Sd = [coef(1) coef(2) coef(3)]
Rd = [coef(4) coef(5) coef(6)]
%
%Espacio de estado observador

%

num = N
roots(N)
den = conv(M,Sd) + conv(N,Rd)

roots(den)
gg = num(4)/den(6)

Sdx = [coef(1) coef(2) coef(3)]
Rdx = [coef(4) coef(5) coef(6)]

Sdxx = 1/2*[coef(1) coef(2) coef(3)]
Rdxx = 1/-30.75*[coef(4) coef(5) coef(6)]

ax1=Rdx(1)
ax2=Rdx(2)
ax3=Rdx(3)

bx1=Sdx(1)
bx2=Sdx(2)
bx3=Sdx(3)


A = [0 0 -bx3;
    1 0 -bx2;
    0 1 -bx1]

B = [ax3 1;
    ax2 0;
    ax1 0]

%%%Lo del doc

numx = conv(N,Rd)
denx = conv(M,Sd)
roots(N)
roots(numx)
roots(M)
roots(denx)
rlocus(numx,denx)
axis([-8 8 -8 8])
axis([-2 2 -2 2])
%ggg= chi(3)/aa;
Ax = [0 -bx3;
    1 -bx2]
poly(Ax)
Sd
pSR = [ax2-ax1*bx2 ax3-ax1*bx3]
pR = (ax1*Sd + conv([0,1],pSR))
Bx = [1 -pSR(2);
    0 -pSR(1)]
Dx = [0 -ax1]
Cx = [0 1]