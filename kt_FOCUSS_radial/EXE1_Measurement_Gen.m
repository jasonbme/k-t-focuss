

%% Add path
clear all;
addpath('bin');
addpath('data');


%% make mask (centered intensity)
DsRate = 6;     % define
disp(['Downsample rate: ',num2str(DsRate)]);

%% Select Reference type for ME/MC
Ref_Type = 'FullSingleFrame';
% Ref_Type = 'DiastoleFrames';

%% Read full measurement

% % Our Simulation Set-up % %
load('full_rad.mat'); % file name => orig (r, theta, t)
Scaling = 1e6;
orig = orig*Scaling;
[DownRadialSino DownRadialMask] = radial_downsampling(orig,DsRate);

% % % Set up your original measurements % % %
% DownRadialSino = (r,phi,frame) = Complex value
% DownRadialMask = (phi, frame) = Degree (0~180)

%% Radial -> Cartesian
[DownCartSino mask] = radial2cart(DownRadialSino, DownRadialMask);

%% you should set diastole frames or a single full frame for making the reference image for ME/MC procedure.

if strcmp(Ref_Type,'DiastoleFrames') == 1
    
    Diastoleframes = [27:32]; % define
    
    sum_down_cart = sum(DownCartSino(:,:,Diastoleframes),3);
    sum_mask = sum(mask(:,:,Diastoleframes),3);
    
    ref_sino = sum_down_cart./sum_mask;
    ref_sino(isnan(ref_sino)) = 0;
    ref_sino(isinf(ref_sino)) = 0;
    
    % % change inverse function by your measurements
    ref = fftshift(fftshift(ifft2(ifftshift(ifftshift(ref_sino,1),2)),2),1);
    
    figure(92);
    imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
    title('Reference image for ME/MC using diastole frames');
    pause(1);
%     close(92);
    
elseif strcmp(Ref_Type,'FullSingleFrame') == 1
    
    Ref_frame = 1; % define
    [nR nA] = size(orig(:,:,Ref_frame));
    Deg = zeros(nA,1);
    Deg(:) = [0:nA-1]/nA*180;
    
    % % full sampling
    [ref_sino full_mask] = radial2cart(orig(:,:,Ref_frame), Deg);
    
    % % change inverse function by your measurements
    ref = fftshift(fftshift(ifft2(ifftshift(ifftshift(ref_sino,1),2)),2),1);
    
    figure(92);
    imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
    title('Reference image for ME/MC using a full single frame');
    pause(1);
    close(92);
end

%% save
save Y.mat DownRadialSino DownRadialMask DownCartSino mask ref
