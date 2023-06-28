clc 
clear all 
close all 

syms phi theta psi

C1 = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
C2 = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
C3 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];

C_BN = C1*C2*C3;
pretty(C_BN)
%%
%____________________________________

C_BN_t = [cos(psi)*cos(theta)  cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi)     sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta); ....
          cos(theta)*sin(psi)  cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)     cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi); ...
          -sin(theta)          cos(theta)*sin(phi)                                  cos(phi)*cos(theta);]
   
r1 = C_BN*C_BN_t
simplify(r1)
