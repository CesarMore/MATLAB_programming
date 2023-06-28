clc
clear all 
close all 

T = 0.01;
t = [0:T:5];
x = [0 0 0];
u = 0;

for i=1:length(t)-1
    [tt,xx] = ode45(@Eje98A,[t(i) t(i+1)],x(i,:),[],u(i));
    x(i+1,:) = xx(end,:);
    
   
   
end