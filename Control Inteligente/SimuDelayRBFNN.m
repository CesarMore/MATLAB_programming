clc
clear all 
close all 
delete('uDNN.txt');
delete('udDNN.txt');
lag = 0.05; %tau
x = [0.1;0.1;0.1;0;0;0;0;0;0];

sol = ode45(@DelayRBFNN,lag,x,[0 20]);







