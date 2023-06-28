syms p r q Ix Iy Iz Ixy Ixz Iyz Izx Izy u v w x y z

S = [0 -r q;
     r 0 -p;
    -q p  0];

Io = [Ix  -Ixy -Ixz;
     -Iyz  Iy  -Iyz;
     -Izx -Izy Iz];
 omega = [p; q; r];
 
 T = S*Io*omega;
 pretty(T)
 
 v = [u;v;w];
 Sr = [0 -z y;
       z 0 -x;
      -y x 0];
  
