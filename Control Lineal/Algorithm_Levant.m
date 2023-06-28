
for k = 1:200
  y(k) = k;     %Trayectoria deseada
end

% codigo

for k = 2:200
    x(k) = x(k-1)+d_y(k-1);       %Integration of X (internal variable) using Euler backward method (x-dot= u)
    u_1(k) = (-miu)*sign(x(k) - y(k-1)); % Formula of u_1 dot in the algorithm 
    U_1(k) = U_1(k-1) + u_1(k);  % Integration of u_1 dot to obtain u_1 using Euler forward method 
    d_y(k) = U_1(k) - lamda*sqrt(abs(x(k) - y(k-1)))*sign(x(k) - y(k-1)); % derivative of y (U in the algorithm)
    
    d_y(k) = y(k) - y(k-1); % Derivative of Y using differece formula
end
