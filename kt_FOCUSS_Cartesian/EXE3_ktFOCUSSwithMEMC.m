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

% % % % % % % % % % % % % % % % % % % % % 
% % Function setting
A = @(x,mask)  fft(x,[],1).*mask;
AT = @(x,mask) ifft(x.*mask,[],1);
ATA = @(x,mask) AT(A(x,mask),mask);

% % % % % % % % % % % % % % % % % % % % % 
% % Parameter setting
Mouter = 2;
Minner = 40;
factor = 0.5;
lambda_focuss = 0;

% % % % % % % % % % % % % % % % % % % % % 
% % init. image
Low_Y = Y;
Low_Y(num_low_phase/2+1:end-num_low_phase/2,:,:) = 0;

%% kt-FOCUSS calculation
X_FOCUSS = KTFOCUSS(A,AT,Y,Low_Y,mask,factor,lambda_focuss, Minner, Mouter);

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
% % Function setting
A = @(x,mask)  fft(x,[],1).*mask;
AT = @(x,mask) ifft(x.*mask,[],1);
ATA = @(x,mask) AT(A(x,mask),mask);

% % % % % % % % % % % % % % % % % % % % % 
% % parameter setting
Mouter = 2;
Minner = 40;
factor = 0.5;
lambda_focuss = 0;
% % % % % % % % % % % % % % % % % % % % % 

Low_Y = Y-A(X_MEMC_tmp,mask);
Low_Y(num_low_phase/2+1:end-num_low_phase/2,:,:) = 0;

X_res = KTFOCUSS(A,AT,Y-A(X_MEMC_tmp,mask),Low_Y,mask,factor,lambda_focuss, Minner, Mouter);

X_MEMC = X_MEMC_tmp+X_res;

%% plot

M = max(abs(X_FOCUSS(:)));

for ff = 1:nframe
    figure(32);
    imagesc(abs(X_FOCUSS(:,:,round(ff))),[0 M]); axis off; axis equal; colormap gray; colorbar; title(['kt-FOCUSS']);
    figure(33);
    imagesc(abs(X_MEMC(:,:,round(ff))),[0 M]); axis off; axis equal; colormap gray; colorbar; title(['kt-FOCUSS with ME/MC']);
    pause(0.01);
end

figure(34);
imagesc(abs(ref),[0 M]); axis off; axis equal; colormap gray; colorbar; title(['Reference image']);

%% save
save Result.mat X_MEMC X_FOCUSS
