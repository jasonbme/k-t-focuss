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

Y = DownCartSino;
mask = double(mask);

[nx ny nframe] = size(Y);

% % % % % % % % % % % % % % % % % % % % % 
% % Function setting 2D-fft
A = @(x,mask) fftshift(fftshift(fft2(x),1),2).*mask;
AT = @(x,mask) ifft2(ifftshift(ifftshift(x.*mask,1),2));
ATA = @(x,mask) AT(A(x,mask2Dt),mask);

% % % % % % % % % % % % % % % % % % % % % 
% % Parameter setting
Mouter = 2;
Minner = 40;
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
for frame=1:nframe
    figure(22);
    imagesc(abs(X_FOCUSS(:,:,frame))); axis off; axis equal; colormap gray; colorbar; title('kt-FOCUSS');
    pause(0.01);
end

close(22);
pause(0.1);

%% ME/MC

disp('ME/MC start: ');

% % % % % % % % % % % % % % % % % % % % % 
% % Parameter setting
px = 3; % patch size of x
py = 3; % patch size of y
ws = 4; % window variation size px+ws*2, ...
% % % % % % % % % % % % % % % % % % % % % 

X_MEMC_tmp = zeros(size(X_FOCUSS));
tic;
for iframe = 1 : nframe
    [motionVect, APRScomputations]=motionEstARPS2(X_FOCUSS(:,:,iframe),ref,px,ws);
    X_MEMC_tmp(:,:,iframe)=motionComp2(ref,motionVect,px);
end

Current_time = toc;
disp(['   ME/MC calculation time - ',num2str(Current_time)]);

%% Residual calculation

% % % % % % % % % % % % % % % % % % % % % 
% % function setting 2D-fft
A = @(x,mask) fftshift(fftshift(fft2(x),1),2).*mask;
AT = @(x,mask) ifft2(ifftshift(ifftshift(x.*mask,1),2));
ATA = @(x,mask) AT(A(x,mask2Dt),mask);

% % % % % % % % % % % % % % % % % % % % % 
% % parameter setting
Mouter = 2;
Minner = 40;
factor = 0.5;
lambda_focuss = 0;
% % % % % % % % % % % % % % % % % % % % % 


X_MEMC_tmp = fftshift(fftshift(X_MEMC_tmp,1),2);

Low_Y = Y-A(X_MEMC_tmp,mask);
LowFreqRatio = 0.05; % change
Low_Y(1:round(nx/2-nx*LowFreqRatio),:,:) = 0;
Low_Y(:,1:round(ny/2-ny*LowFreqRatio),:) = 0;
Low_Y(round(nx/2+nx*LowFreqRatio):end,:,:) = 0;
Low_Y(:,round(ny/2+ny*LowFreqRatio):end,:) = 0;

% X_res = KTFOCUSS(A,AT,Y-A(X_MEMC,mask),Low_Y,mask,factor,lambda_focuss, Minner, Mouter);
X_res = KTFOCUSS_radial(A,AT,Y-A(X_MEMC_tmp,mask),Low_Y,mask,DownRadialSino-cart2radial(A(X_MEMC_tmp,mask),DownRadialMask),DownRadialMask,factor,lambda_focuss, Minner, Mouter);

X_MEMC = fftshift(fftshift(X_MEMC_tmp+X_res,1),2);

%% plot
for ff = 1:nframe
    figure(32);
    imagesc(abs(X_FOCUSS(:,:,round(ff)))); axis off; axis equal; colormap gray; colorbar; title(['kt-FOCUSS']);
    figure(33);
    imagesc(abs(X_MEMC(:,:,round(ff)))); axis off; axis equal; colormap gray; colorbar; title(['kt-FOCUSS with ME/MC']);
    pause(0.01);
end

figure(34);
imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar; title(['Reference image']);

%% save
save Result.mat X_MEMC X_FOCUSS
