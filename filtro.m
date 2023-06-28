clc 
clear all
close all

T = readmatrix("Datos3.txt")
%T1 = readmatrix("Datos3.txt")
figure(1)
X = plot(T(:,1));
hold on 
X = plot(T(:,2))
hold on 
X = plot(T(:,3))
% figure(2)
% X = plot(T(:,2));
% hold on 
% X = plot(-T(:,9))
% figure(3)
% X = plot(T(:,3));
% hold on 
% X = plot(-T(:,10))
% figure(4)
% X = plot(T(:,4));
% figure(5)
% X = plot(T(:,5));
% figure(6)
% X = plot(T(:,6));
% figure(7)
% X = plot(T(:,7));
% figure(8)
% X = plot(T(:,8));
% figure(9)
% X = plot(T(:,9));
% figure(10)
% X = plot(T(:,10));
% figure(11)
% X = plot(T(:,11));
% hold on
% Y = plot(T(:,9));
% 
% figure(2)
% X = plot(T(:,2));
% hold on
% Y = plot(-T(:,8));