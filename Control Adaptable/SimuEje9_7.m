clc 
%clear all 
%close all

T = 0.01;
t = [0:T:10];

x = 5;
u1 = 0;
u2 = 0;
alpha = sqrt(2);

for i=1:length(t)-1
    [tt,xx] = ode45(@Eje9_7,[t(i) t(i+1)],x(i,:),[],u1(i));
    x(i+1,:) = xx(end,:);
    
    u1(i+1) = -x(i+1,1) + (x(i+1,1))^3 -alpha*x(i+1,1);
    u2(i+1) = -x(i+1,1) + (x(i+1,1))^3 - x(i+1,1)*sqrt( (1-(x(i+1,1))^2)^2 +1);
    
end
% 
% figure(1)
% plot(t,x)
% hold on
% plot(t,u1)
% legend('x','u1')
% xlabel('tiempo [s]')
% ylabel('x')
% 
% figure(2)
% plot(t,x)
% hold on
% plot(t,u2)
% legend('x','u2')
% xlabel('tiempo [s]')
% ylabel('x')
% ylim([-2 2])
% grid on 

figure(1)
plot(out.tout,out.x)
grid on
legend('x')
figure(2)
plot(out.tout,out.u1)
grid on
legend('u1')
