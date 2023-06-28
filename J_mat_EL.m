clear all
clc
syms Ixx Iyy Izz theta psi phi i a s p  x y z

W=[1        0            -sin(theta);
    0    cos(phi)    sin(phi)*cos(theta);
    0   -sin(phi)   cos(phi)*cos(theta)]

I=[Ixx 0 0;
    0 Iyy 0;
    0 0 Izz]

J=transpose(W)*I*W

simplify(J)


eta = [ phi*(p), theta*(p), psi*(p)]
T= eta*J*transpose(eta)

%%
%%%%%%%%%%%
clear all 
clc
syms x y z a
A=[1 3 -2 x;
   0 5 -1 y;
   0 15 -3 z]

R=rref(A)

B=[0 -5 1 2;
   3 -1 1 3]
Rr=rref(B)

C=[3 -1 a ;
   2 -1 5 ;
   1 -1 2 ]
rr=rref(C)

M=[-1/15 1/3;
    -1/5 0]
m=[1 2;
    1 3]
M*m