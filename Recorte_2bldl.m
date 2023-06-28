%con el cursor 
close all 
clear all
clc
i = imread('Documents/MATLAB/imagesP1/bld1.jpg');
figure(1)
imshow(i)
[rec, pos] = imcrop(i)
figure(2)
imshow(rec)


