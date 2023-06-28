clc
clear all
close all

ti = 0;
tf = 5;
T =0.01;
t = [ti:T:tf];
c = [-0.8 -0.4 0 0.4 0.8]; %centro
b = 0.25; %bias
x = [1 0.3];
theta = rand(5,1); %pesos iniciales
psi = 0.5; %cota de la perturbacion 
u = 0;
%zeta = 0;
epsilon = 0.1;

Gamma = diag([0.1 0.1 0.1 0.1 0.1]);
sigma = 0.1;
sigma_psi = 0.1;
gamma = 0.1;

for i=1:length(t)
    [xx,tt] = ode45(@modelo1, [t(i) t(i+1)], x(:,1),[],u(i));
    x(i+1,:) = xx(end,:);
    
    x1 = x(i+1,1);

    for j = 1:5
        zeta(j,1) = exp(abs(x1 - c(j))^2 / b);
    end
    beta1 = psi *tanh(x1/epsilon);
    alpha = -x(i+1,1) - theta'*zeta - beta1; 
    
    PARCIAL = 1;

    z1 = x1;
    z2 = x(i+1,2) - alpha;

    %u
    u = -z1 - z2 -1 -theta  ;

    thetap = theta + T*(Gamma*(zeta*(z1-z2*PARCIAl)-sigma*theta));
    psip = psi + T*(gamma*(z1*tanh(z1/epsilon) + z2*PARCIAL*tanh(z2*PARCIAL/epsilon) - sigma_psi*(psi)));

end