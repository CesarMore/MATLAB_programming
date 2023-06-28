clc
clear all
close all
warning off
%k=imread('C:\Users\cmdjd\OneDrive - CINVESTAV\Documentos\MATLAB\imagesP1\wdg3.jpg'); %cargar imagen 
k=imread('C:\Users\cmdjd\OneDrive - CINVESTAV\Documentos\MATLAB\imagesP1\step2line.jpg'); % cargar image

figure;
imshow(k), title('Imagen original');

x=rand(size(k));
k(x(:)>0.95)=255;

sto=[];
[a b]=size(k);
output=zeros(a,b);
for i=2:a-1
    for j=2:b-1
        sto=[k(i-1,j-1),k(i-1,j),k(i-1,j+1),k(i,j-1),k(i,j)...
            ,k(i,j+1),k(i+1,j-1),k(i+1,j),k(i+1,j+1)];
        es=min(sto);
            output(i,j)=es;
        sto=[];
    end
end
% figure;
% imshow(uint8(output)), title('');



x=rand(size(k));
k(x(:)<0.05)=0;

sto=[];
[a b]=size(k);
output=zeros(a,b);
for i=2:a-1
    for j=2:b-1
        sto=[k(i-1,j-1),k(i-1,j),k(i-1,j+1),k(i,j-1),k(i,j)...
            ,k(i,j+1),k(i+1,j-1),k(i+1,j),k(i+1,j+1)];
        es=max(sto);
            output(i,j)=es;
        sto=[];
    end
end
figure;
imshow(uint8(output)), title('Imagen binarizada');