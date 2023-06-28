clc, clear;

ti = 0;
tf = 80;
Ts = 0.01;
t = ti:Ts:tf;

% condiciones iniciales
q   = [-0.5, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
Gamma = [0, 0, 0, 0];

% parametros del sistema
m   = 0.460;
g   = 9.81;
Ixx = 2.24e-3;
Iyy = 2.90e-3;
Izz = 5.30e-3;

%% referencias

% z deseada
zd   = 2*ones(1,length(t));
zdp  = zeros(1,length(t));
zdpp = zeros(1,length(t));

% psi deseada
psid   = zeros(1,length(t));
psidp  = zeros(1,length(t));
psidpp = zeros(1,length(t));

% x-y deseadas
t1 = t(1:2000);
t2 = t(1:2001);

% coordenadas x
[xd1,xdp1,xdpp1] = segmento(t1,0,-10);
xd(1:2000) = xd1;       xdp(1:2000) = xdp1;     xdpp(1:2000) = xdpp1;    
xd(2001:4000) = -10;    xdp(2001:4000) = 0;     xdpp(2001:4000) = 0;
[xd3,xdp3,xdpp3] = segmento(t1,-10,10);
xd(4001:6000) = xd3;    xdp(4001:6000) = xdp3;  xdpp(4001:6000) = xdpp3;
xd(6001:8001) = 0;      xdp(6001:8001) = 0;     xdpp(6001:8001) = 0;

% coordenadas y
yd(1:2000) = 0;         ydp(1:2000) = 0;        ydpp(1:2000) = 0;
[yd2,ydp2,ydpp2] = segmento(t1,0,10);
yd(2001:4000) = yd2;    ydp(2001:4000) = ydp2;  ydpp(2001:4000) = ydpp2;  
yd(4001:6000) = 10;     ydp(4001:6000) = 0;     ydpp(4001:6000) = 0;
[yd4,ydp4,ydpp4] = segmento(t2,10,-10);
yd(6001:8001) = yd4;    ydp(6001:8001) = ydp4;  ydpp(6001:8001) = ydpp4;


%% Ganancias de control

k1_x = 10;     k1_y = 10;       k1_z = 20;
k2_x = 3;      k2_y = 3;        k2_z = 7;

k1_phi = 50;   k1_theta = 50;   k1_psi = 50;
k2_phi = 10;   k2_theta = 10;   k2_psi = 10;


%% Integracion numerica
for i=1:length(t)-1
    
    [tt,qq] = ode45(@mod_cuadri, [t(i), t(i+1)], q(i,:), [], Gamma(i,:));
    q(i+1,:) = qq(end,:);
    
    % variables de estado
    x  = q(i+1,1);  y  = q(i+1,2);  z  = q(i+1,3);  phi  = q(i+1,4);   theta  = q(i+1,5);   psi  = q(i+1,6);    
    xp = q(i+1,7);  yp = q(i+1,8);  zp = q(i+1,9);  phip = q(i+1,10);  thetap = q(i+1,11);  psip = q(i+1,12);
    
    
    % control de altitud --------------------------------------------------
    ez  = zd(i+1)  - z; 
    ezp = zdp(i+1) - zp;
    
    % F
    Gamma(i+1,1) = (m/(cos(phi)*cos(theta)))*(zdpp(i+1) + k2_z*ezp + k1_z*ez + g);
    
    
    % control de psi ------------------------------------------------------
    epsi =  psid(i+1)  - psi;
    epsip = psidp(i+1) - psip;
    
    % tau_psi
    Gamma(i+1,4) = Izz*(psidpp(i+1) + k2_psi*epsip + k1_psi*epsi);
    
    
    % control de x-y ------------------------------------------------------
    ex  = xd(i+1)  - x; 
    exp = xdp(i+1) - xp;
    
    ey  = yd(i+1)  - y; 
    eyp = ydp(i+1) - yp;
    
    mu_x = (m/Gamma(i+1,1))*(xdpp(i+1) + k2_x*exp + k1_x*ex);
    mu_y = (m/Gamma(i+1,1))*(ydpp(i+1) + k2_y*eyp + k1_y*ey);
    
    var_x(i+1) =  mu_x*sin(psi) - mu_y*cos(psi);
    var_y(i+1) = (mu_x*cos(psi) + mu_y*sin(psi))/cos(phi);
    
    phid(i+1)   = asin(var_x(i+1));
    thetad(i+1) = asin(var_y(i+1));
    
    phidp(i+1) = (phid(i+1) - phid(i))/Ts;
    thetadp(i+1) = (thetad(i+1) - thetad(i))/Ts;
    
    phidpp(i+1) = (phidp(i+1) - phidp(i))/Ts;
    thetadpp(i+1) = (thetadp(i+1) - thetadp(i))/Ts;
    
    
    ephi  = phid(i+1)  - phi;
    ephip = phidp(i+1) - phip;
    
    etheta  = thetad(i+1)  - theta;
    ethetap = thetadp(i+1) - thetap;
     
    % tau_phi
    Gamma(i+1,2) = Ixx*(phidpp(i+1) + k2_phi*ephip     + k1_phi*ephi);
    % tau_theta
    Gamma(i+1,3) = Iyy*(thetadpp(i+1) + k2_theta*ethetap + k1_theta*etheta);

end


%% Plotting

% Estados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure(1);
f1.Position = [680, 558, 1200, 550];
subplot(2,2,1)
plot(t, xd, 'k--', t, q(:,1), 'r', t, yd, 'k--', t, q(:,2), 'b','LineWidth', 2), grid on
set(gca, 'FontSize', 12)
legend('$x_d$', '$x$', '$y_d$', '$y$','Interpreter', 'latex', 'Location', 'SouthEast')
xlabel('[seg]')
ylabel('[m]')

subplot(2,2,3)
plot(t, xdp, 'k--', t, q(:,7), 'r', t, ydp, 'k--', t, q(:,8), 'b', 'LineWidth', 2), grid on
set(gca, 'FontSize', 12)
legend('$\dot{x}_d$', '$\dot{x}$', '$\dot{y}_d$', '$\dot{y}$','Interpreter', 'latex', 'Location', 'NorthEast')
xlabel('[seg]')
ylabel('[m/s]')


subplot(2,2,[2,4])
plot3(xd, yd, zd, 'k--', q(:,1), q(:,2), q(:,3), 'g','LineWidth', 2), grid on
set(gca, 'FontSize', 12)
legend('ref', 'cuadri', 'Location', 'NorthEast')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xlim([-12, 2]);
ylim([-2, 12]);
zlim([0, 2.2]);

% Estados orientacion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure(2);
f2.Position = [680, 558, 600, 250];
plot(t, q(:,4), 'r', t, q(:,5), 'b', t, q(:,6), 'm', 'LineWidth', 2), grid on
set(gca, 'FontSize', 12)
legend('\phi', '\theta', '\psi')
xlabel('[seg]')
ylabel('[rad]')
%ylim([-0.2, 0.2]);


% Entradas de control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure(3);
f3.Position = [680, 558, 600, 350];
subplot(2,1,1)
plot(t,Gamma(:,1), 'r', 'LineWidth', 2), grid on
set(gca, 'FontSize', 12)
legend('f_q')
xlabel('[seg]')
ylabel('[N]')

subplot(2,1,2)
plot(t,Gamma(:,2), 'r', t,Gamma(:,3), 'b', t,Gamma(:,4), 'm', 'LineWidth', 2), grid on
set(gca, 'FontSize', 12)
legend('\tau_{\phi}', '\tau_{\theta}', '\tau_{\psi}')
xlabel('[seg]')
ylabel('[Nm]')
ylim([-0.005, 0.005]);

