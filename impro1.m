clc
clear all
close all
I = imread('rice.jpg');
J = imread('cameraman.tif');
for i = 1:256
    for j = 1:256
        output(i,j)=(I(i,j)+J(i,j)/2);
    end
end
imshow(output);
