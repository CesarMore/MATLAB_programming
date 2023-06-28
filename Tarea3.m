clear all 
close all
clc
colr_1 = [104 201 80];
colr_2 = [95 64 217];
img = imread("Lena.png");

l1 = round(linspace(colr_1(1),colr_2(1),512));
l2 = round(linspace(colr_1(2),colr_2(2),512));
l3 = round(linspace(colr_1(3),colr_2(3),512));

img = uint8(img);
for x=1:length(l1)
    img_1(:,x,1)=l1(x);
    img_2(:,x,2)=l2(x);
    img_3(:,x,3)=l3(x);
    figure(1),imshow(img_1)
    drawnow
end

    
    
