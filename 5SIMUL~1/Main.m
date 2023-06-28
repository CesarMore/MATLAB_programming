clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIEMPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf = 40;             % Tiempo de simulacion en segundos (s)
ts = 0.1;            % Tiempo de muestreo en segundos (s)
t = 0: ts: tf;       % Vector de tiempo
N = length(t);       % Muestras
%%%%%%%%%%%%%%%%%%%%%%%% PARAMETROS DEL ROBOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0.2; % Distancia hacia el punto de control en metros (m)

%%%%%%%%%%%%%%%%%%%%%%%% CONDICIONES INICIALES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = zeros (1,N+1); % Posición en el centro del eje que une las ruedas (eje x) en metros (m)
y1 = zeros (1,N+1); % Posición en el centro del eje que une las ruedas (eje y) en metros (m)

x1(1) = 1;          % Posicion inicial eje x
y1(1) = -1;        % Posicion inicial eje y
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODELO GEOMETRICO %%%%%%%%%%%%%%%%%%%%%%%%%%%%

hx = zeros(1, N+1);  % Posicion en el punto de control (eje x) en metros (m)
hy = zeros(1, N+1);  % Posicion en el punto de control (eje y) en metros (m)
phi = zeros(1, N+1); % Orientacion del robot en radianes (rad)

phi(1) = 0;   % Orientacion inicial del robot

hx(1) = x1(1) + a*cos(phi(1)); % Posicion en el punto de control del robot en el eje x
hy(1) = y1(1) + a*sin(phi(1)); % Posicion en el punto de control del robot en el eje y


%%%%%%%%%%%%%%%%%%%%%% VELOCIDADES DE REFERENCIA %%%%%%%%%%%%%%%%%%%%%%%%%%

uRef = 0.2*ones(1,N); 
wRef = 0.1*ones(1,N); 


%%%%%%%%%%%%%%%%%%%%%%%%% BUCLE DE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:N 
    
    %%%%%%%%%%%%%%%%%%%%% MODELO CINEMATICO %%%%%%%%%%%%%%%%%%%%%%%%%

    phi(k+1)=phi(k)+wRef(k)*ts;

    xp1=uRef(k)*cos(phi(k+1));
    yp1=uRef(k)*sin(phi(k+1));

    x1(k+1)=ts*xp1+ x1(k);
    y1(k+1)=ts*yp1+ y1(k);

    hx(k+1)=x1(k+1)+a*cos(phi(k+1)); 
    hy(k+1)=y1(k+1)+a*sin(phi(k+1));
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULACION VIRTUAL 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a) Configuracion de escena

scene=figure;  % Crear figura (Escena)
set(scene,'Color','white'); % Color del fondo de la escena
set(gca,'FontWeight','bold') ;% Negrilla en los ejes y etiquetas
sizeScreen=get(0,'ScreenSize'); % Retorna el tamaño de la pantalla del computador
set(scene,'position',sizeScreen); % Congigurar tamaño de la figura
axis equal; % Establece la relación de aspecto para que las unidades de datos sean las mismas en todas las direcciones.
axis([-5 5 -5 5 0 1]); % Ingresar limites minimos y maximos en los ejes x y z [minX maxX minY maxY minZ maxZ]
view([-40 30]); % Orientacion de la figura
grid on; % Mostrar líneas de cuadrícula en los ejes
box on; % Mostrar contorno de ejes
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)'); % Etiqueta de los ejes
camlight right % Luz para la escena

% b) Graficar robots en la posicion inicial
MobileRobot;
H1=MobilePlot(x1(1),y1(1),phi(1));hold on;

% c) Graficar Trayectorias
H2=plot3(hx(1),hy(1),0,'r','lineWidth',2);

% d) Bucle de simulacion de movimiento del robot

step=10; % pasos para simulacion

for k=1:step:N

    delete(H1);    
    delete(H2);
    
    H1=MobilePlot(x1(k),y1(k),phi(k));
    H2=plot3(hx(1:k),hy(1:k),zeros(1,k),'r','lineWidth',2);
    
    pause(ts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graficas %%%%%%%%%%%%%%%%%%%%%%%%%%%%
graph=figure;  % Crear figura (Escena)
set(graph,'position',sizeScreen); % Congigurar tamaño de la figura
subplot(211)
plot(t,uRef,'b','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('m/s');
subplot(212)
plot(t,wRef,'r','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('[rad/s]');


