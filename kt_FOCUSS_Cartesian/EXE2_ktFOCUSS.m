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

Y = DownSino;
mask = double(mask);

[nx ny nframe] = size(Y);

% % % % % % % % % % % % % % % % % % % % change
% % function setting
A = @(x,mask)  fft(x,[],1).*mask;
AT = @(x,mask) ifft(x.*mask,[],1);
ATA = @(x,mask) AT(A(x,mask),mask);

% % % % % % % % % % % % % % % % % % % % 
% % parameter setting
Mouter = 2;
Minner = 40;
factor = 0.5;
lambda_focuss = 0;
% % % % % % % % % % % % % % % % % % % %

Low_Y = Y;
Low_Y(num_low_phase/2+1:end-num_low_phase/2,:,:) = 0;

X_FOCUSS = KTFOCUSS(A,AT,Y,Low_Y,mask,factor,lambda_focuss, Minner, Mouter);

%% plot
for frame=1:nframe
    figure(22);
    imagesc(abs(X_FOCUSS(:,:,frame))); axis off; axis equal; colormap gray; colorbar; title('kt-FOCUSS');
    pause(0.01);
end
