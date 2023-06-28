   
for k=1:200
  f(k)=k*sin(k);     % Trayectoria deseada
end

U_1 = 1;
u_1 = 1;
    % codigo
for k=2:200
    U_1(k) = U_1(k-1)+u_1(k);  % Integration of u_1 dot to obtain u_1 using Euler forward method 
    d_U(k) = U_1(k) - lamda*sqrt(abs(x(k)-f(k-1)))*sign(x(k)-f(k-1)); % derivative of U in the algorithm)
    x(k) = x(k-1) + d_f(k-1);       %Integration of X (internal variable) using Euler backward method (x-dot= u)
    u_1(k) = (-miu)*sign(x(k)-f(k-1)); % Formula of u_1 dot in the algorithm 

    d_f(k)=f(k)-f(k-1); % Derivative of f using differece formula
end
