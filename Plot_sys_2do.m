clear all 
close all
clc

x0 = [5; 5]; % condiciones iniciales
t = [0:0.05:60]; % tiempo

for i=1:length(t)-1            
    [tt, xx] = ode45(@my_system, [t(i) t(i+1)], x0);
    
    x(i+1,:) = xx(end,:);
    
    x1real = x(i+1,1);
    x2real = x(i+1,2);
    x1p = sin(x1real - x2real) + x1real^2 + x2real^2*sin(t);
    x2p = log(1 + x1real + x2real) - x2real^2;
end


figure(1);
plot(t, x(:,1), 'b');
legend('x1');
xlabel('Time');
ylabel('State Variable');
title('Simulation');
figure(2);
plot(t, x(:,2), 'r');
legend('x2');
xlabel('Time');
ylabel('State Variable');
title('Simulation');
%x punto
% figure(3);
% plot(t, x1p(i,:), 'b');
% legend('x1');
% xlabel('Time');
% ylabel('State Variable');
% title('Simulation');
% figure(4);
% plot(t, x2p(i,:), 'r');
% legend('x2');
% xlabel('Time');
% ylabel('State Variable');
% title('Simulation');
% x1p vs x2p
figure(5);
plot(x1p,'b', x2p, 'r');
legend('x2');
xlabel('Time');
ylabel('State Variable');
title('Simulation');


function dxdt = my_system(t, x)
    x1 = x(1);
    x2 = x(2);
    
    %dxdt = zeros(2,1);
    dxdt(1,1) = sin(x1 - x2) + x1^2 + x2^2*sin(t);
    dxdt(2,1) = log(1 + x1 + x2) - x2^2;
end
