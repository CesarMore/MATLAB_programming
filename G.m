function G
v1 = input('Dame el vector v1');
v2 = input('Dame el vector v2');
v3 = input('Dame el vector v3');


A=[v1;v2;v3];

while(det(A)==0)
disp('Estos vectores no son linealmente independientes');
    
v1 = input('Dame el vector v1');
v2 = input('Dame el vector v2');
v3 = input('Dame el vector v3');
    

A=[v1;v2;v3];
end

u1=v1;
u2= v2-(dot(v2,u1)/norm(u1)^2)*u1;
u3= v3-(dot(v3,u2)/norm(u2)^2)*u2-(dot(v3,u1)/norm(u1)^2)*u1;


disp(strcat('El vector u1 es', mat2str(u1)));

disp(strcat('El vector u2 es', mat2str(u2)));

disp(strcat('El vector u3 es', mat2str(u3)));

end