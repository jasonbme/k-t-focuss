
%% Clear all
clear all
close all
clc;

addpath('bin');
addpath('data');

%% Read Measurement
load Y.mat

%% kt-FOCUSS

disp('kt-FOCUSS start: ');

% Y = DownRadialSino;
% mask = double(DownRadialMask);

[Y mask] = radial2cart(DownRadialSino, DownRadialMask);

[nx ny nframe] = size(Y);

% % % % % % % % % % % % % % % % % % % % % 
% % Function setting 2D-fft
A = @(x,mask) fftshift(fftshift(fft2(x),1),2).*mask;
AT = @(x,mask) ifft2(ifftshift(ifftshift(x.*mask,1),2));
ATA = @(x,mask) AT(A(x,mask2Dt),mask);

% % % % % % % % % % % % % % % % % % % % % 
% % Parameter setting
Mouter = 2;
Minner = 20;
factor = 0.5;
lambda_focuss = 0;

% % init. image
Low_Y = Y;
LowFreqRatio = 0.05; % change
Low_Y(1:round(nx/2-nx*LowFreqRatio),:,:) = 0;
Low_Y(:,1:round(ny/2-ny*LowFreqRatio),:) = 0;
Low_Y(round(nx/2+nx*LowFreqRatio):end,:,:) = 0;
Low_Y(:,round(ny/2+ny*LowFreqRatio):end,:) = 0;


%% kt-FOCUSS calculation

X_FOCUSS = KTFOCUSS_radial(A,AT,Y,Low_Y,mask,DownRadialSino,DownRadialMask,factor,lambda_focuss, Minner, Mouter);
X_FOCUSS = fftshift(fftshift(X_FOCUSS,1),2);


%% plot
M = max(abs(X_FOCUSS(:)));
for frame=1:nframe
    figure(22); imagesc(abs(X_FOCUSS(:,:,frame)),[0 M]); axis off; axis equal; colormap gray; colorbar; title('kt-FOCUSS');
    pause(0.001);
end
