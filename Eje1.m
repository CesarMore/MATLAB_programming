clc 
clear
 
syms t s
 
a=1
b=1
A = [0 1 0 0 0;-a^2 0 0 0 0;0 0 0 1 0;0 0 0 0 1;0 0 -a^3 -3*a^2 -3*a]
bv = [0;0;0;0;1]
c = [1;0;-b^2;0;1]

M= [s*eye(5)- A]
 
pi(s) = det(M)
 
vp = eig(A) 
%%
sigma = [M, bv; -transpose(c), 0]

det1 = det(sigma)
ceros = roots([1 0 0 0 -1])

