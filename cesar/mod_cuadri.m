function dqdt = mod_cuadri(t, q, Gamma)   
    % Variables de estado    
    x  = q(1);  y  = q(2);  z  = q(3);  phi  = q(4);   theta  = q(5);   psi  = q(6);    
    xp = q(7);  yp = q(8);  zp = q(9);  phip = q(10);  thetap = q(11);  psip = q(12);
    
    xi   = [x; y; z];
    eta  = [phi; theta; psi]; 
    xip  = [xp; yp; zp];
    etap = [phip; thetap; psip]; 
    
    % Entradas de control
    F         = Gamma(1);
    tau_phi   = Gamma(2);
    tau_theta = Gamma(3);
    tau_psi   = Gamma(4);

    % parametros del sistema
    m   = 0.460;
    g   = 9.81;
    Ixx = 2.24e-3;
    Iyy = 2.90e-3;
    Izz = 5.30e-3;
    
    %% Modelo dinamico
    
    % Traslacion ----------------------------------------------------------
    M_xi = m*eye(3);
    
    G_xi = [0;0;m*g];
      
    % Rotacion ------------------------------------------------------------
    M11_eta = Ixx;
    M12_eta = 0;
    M13_eta = -Ixx*sin(theta);
    M21_eta = 0;
    M22_eta = Iyy*cos(phi)^2 + Izz*sin(phi)^2;
    M23_eta = (Iyy-Izz)*sin(phi)*cos(phi)*cos(theta);
    M31_eta = -Ixx*sin(theta);
    M32_eta = (Iyy-Izz)*sin(phi)*cos(phi)*cos(theta);
    M33_eta = Ixx*sin(theta)^2 + Iyy*sin(phi)^2*cos(theta)^2 + Izz*cos(phi)^2*cos(theta)^2;
    
    M_eta = [M11_eta, M12_eta, M13_eta;
             M21_eta, M22_eta, M23_eta;
             M31_eta, M32_eta, M33_eta];

    C11_eta = 0;
    C12_eta = (Iyy-Izz)*(psip*sin(phi)^2*cos(theta) - psip*cos(phi)^2*cos(theta) + thetap*sin(phi)*cos(phi)) - Ixx*psip*cos(theta);
    C13_eta = (Izz-Iyy)*psip*sin(phi)*cos(phi)*cos(theta)^2;
    C21_eta = (Izz-Iyy)*(psip*sin(phi)^2*cos(theta) - psip*cos(phi)^2*cos(theta) + thetap*sin(phi)*cos(phi)) + Ixx*psip*cos(theta);
    C22_eta = (Izz-Iyy)*phip*sin(phi)*cos(phi);
    C23_eta = -(Ixx - Iyy*sin(phi)^2 - Izz*cos(phi)^2)*psip*sin(theta)*cos(theta);
    C31_eta = (Iyy-Izz)*psip*sin(phi)*cos(phi)*cos(theta)^2 - Ixx*thetap*cos(theta);
    C32_eta = (Izz-Iyy)*(phip*sin(phi)^2*cos(theta) - phip*cos(phi)^2*cos(theta) + thetap*sin(phi)*cos(phi)*sin(theta))... 
        + (Ixx - Iyy*sin(phi)^2 - Izz*cos(phi)^2)*psip*sin(theta)*cos(theta);
    C33_eta = (Iyy-Izz)*phip*sin(phi)*cos(phi)*cos(theta)^2 + (Ixx - Iyy*sin(phi)^2 - Izz*cos(phi)^2)*thetap*sin(theta)*cos(theta);

    C_eta = [C11_eta, C12_eta, C13_eta;
             C21_eta, C22_eta, C23_eta;
             C31_eta, C32_eta, C33_eta];

     
    % Fuerzas externas aplicadas ------------------------------------------
    fI = [cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
          cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
          cos(phi)*cos(theta)]*F;
    
    tau = [tau_phi;
           tau_theta;
           tau_psi];
       
    %% Ecuaciones de estado
    dqdt(1:3,1)   = xip;
    dqdt(4:6,1)   = etap;
    dqdt(7:9,1)   = inv(M_xi)*(fI - G_xi);
    dqdt(10:12,1) = inv(M_eta)*(tau - C_eta*etap);

end