clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIEMPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf = 40;             % Tiempo de simulacion en segundos (s)
ts = 0.1;            % Tiempo de muestreo en segundos (s)
t = 0: ts: tf;       % Vector de tiempo
N = length(t);       % Muestras

%%%%%%%%%%%%%%%%%%%%%%%% PARAMETROS DEL ROBOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.1;

%%%%%%%%%%%%%%%%%%%%%%%% CONDICIONES INICIALES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = zeros(1, N+1);  % Posicion en el centro del robot (eje x) en metros (m)
y1 = zeros(1, N+1);  % Posicion en el centro del robot (eje y) en metros (m)
z1 = zeros(1, N+1);  % Posicion en el centro del robot (eje z) en metros (m)

phi = zeros(1, N+1); % Orientacion del robot en radianes (rad)

x1(1)=0; % Posicion inicial eje x
y1(1)=0; % Posicion inicial eje y
z1(1)=0; % Posicion inicial eje y

phi(1)=0; % Orientacion inicial del robot

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUNTO DE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%

hx = zeros(1, N+1);  % Posicion en el punto de control (eje x) en metros (m)
hy = zeros(1, N+1);  % Posicion en el punto de control (eje y) en metros (m)
hz = zeros(1, N+1);  % Posicion en el punto de control (eje z) en metros (m)

hx(1) = x1(1)+a*sin(phi(1)); % Posicion en el punto de control del robot en el eje x
hy(1) = y1(1)-a*cos(phi(1)); % Posicion en el punto de control del robot en el eje y
hz(1) = z1(1); % Posicion en el punto de control del robot en el eje z


%%%%%%%%%%%%%%%%%%%%%% VELOCIDADES DE REFERENCIA %%%%%%%%%%%%%%%%%%%%%%%%%%

uf =  0.1*ones(1,length(t)); % velocidad lineal frontal (eje x)
ul =  0.2*ones(1,length(t)); % velocidad lineal lateral (eje y)
uz =  0.08*ones(1,length(t));  % velocidad lineal en z 
w  = -0.01*ones(1,length(t)); % velocidad angular


%%%%%%%%%%%%%%%%%%%%%%%%% BUCLE DE SIMULACION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(t)

    %%%%%%%%%%%%%%%%%%%%% MODELO CINEMATICO %%%%%%%%%%%%%%%%%%%%%%%%%
    x1p = uf(k)*cos(phi(k))-ul(k)*sin(phi(k));
    y1p = uf(k)*sin(phi(k))+ul(k)*cos(phi(k));
    z1p = uz(k);
    phip = w(k);
    
    % Integral numérica (método de Euler)
    x1(k+1)=x1(k)+ts*x1p;
    y1(k+1)=y1(k)+ts*y1p;
    z1(k+1)=z1(k)+ts*z1p;
    phi(k+1)=phi(k)+phip;
    

    hx(k+1) = x1(k+1)+a*sin(phi(k+1)); % Posicion en el punto de control del robot en el eje x
    hy(k+1) = y1(k+1)-a*cos(phi(k+1)); % Posicion en el punto de control del robot en el eje y
    hz(k+1) = z1(k+1); % Posicion en el punto de control del robot en el eje z

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULACION VIRTUAL 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a) Configuracion de escena

scene=figure;  % Crear figura (Escena) 
set(scene,'Color','white'); % Color del fondo de la escena
set(gca,'FontWeight','bold') ;% Negrilla en los ejes y etiquetas
sizeScreen=get(0,'ScreenSize'); % Retorna el tamaño de la pantalla del computador
set(scene,'position',sizeScreen); % Congigurar tamaño de la figura
axis equal; % Establece la relación de aspecto para que las unidades de datos sean las mismas en todas las direcciones.
grid on; % Mostrar líneas de cuadrícula en los ejes
box on; % Mostrar contorno de ejes
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)'); % Etiqueta de los eje

view([-40 30]); % Orientacion de la figura
axis([-5 5 -5 5 0 5]); % Ingresar limites minimos y maximos en los ejes x y z [minX maxX minY maxY minZ maxZ]

% b) Graficar robots en la posicion inicial
H1=Plot_Drone(x1(1),y1(1),z1(1),0,0,phi(1),1); hold on;


% c) Graficar Trayectorias
H2=plot3(hx(1),hy(1),hz(1),'b','LineWidth',2); 


% d) Bucle de simulacion de movimiento del robot

step=10; % pasos para simulacion

for k=1:step:N
    
    delete (H1)
    delete (H2)
    H1=Plot_Drone(x1(k),y1(k),z1(k),0,0,phi(k),1);
    H2=plot3(hx(1:k),hy(1:k),hz(1:k),'b','LineWidth',2);
    
    pause(ts)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graficas %%%%%%%%%%%%%%%%%%%%%%%%%%%%
graph=figure;  % Crear figura (Escena)
set(graph,'position',sizeScreen); % Congigurar tamaño de la figura
subplot(211)
plot(t,uf,'b','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('[m/s]'),hold on;
plot(t,ul,'y','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('[m/s]');
plot(t,uz,'r','LineWidth',2),grid('on'),legend('uf','ul','uz');
subplot(212)
plot(t,w,'r','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('[rad/s]'),legend('w');




