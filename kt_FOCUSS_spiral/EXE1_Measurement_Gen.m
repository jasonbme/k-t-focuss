%% Add path
clear all;
addpath('bin');
addpath('data');


%% make mask (centered intensity)
DsRate = 4;     % define
disp(['Downsample rate: ',num2str(DsRate)]);

%% Select Reference type for ME/MC
Ref_Type = 'FullSingleFrame';
% Ref_Type = 'SumSeveralFrames';

%% Read full measurement
load('full_spiral.mat'); % file name => orig (data, line, t); pos => real = x, imag = y, t
Scaling = 1;
orig = orig*Scaling;

% % define cartesian size
nx = 192;
ny = 192;

[ndata nline nframe] = size(orig);

Data = zeros(ndata, nline/DsRate, nframe);
Pos = zeros(ndata, nline/DsRate, nframe);

% % downsampling
for iframe = 1: nframe
    Data(:,:,iframe) = orig(:,1+mod(iframe-1,DsRate):DsRate:nline,iframe);
    Pos(:,:,iframe) = pos(:,1+mod(iframe-1,DsRate):DsRate:nline,iframe);
end

[DownSpiralSino DownMask] = spiral2cart(Data,Pos,nx,ny); 
%% you should set diastole frames or a single full frame for making the reference image for ME/MC procedure.

if strcmp(Ref_Type,'SumSeveralFrames') == 1
    
    SumSeveralFrames = [20:40]; % define
    
    sum_down_cart = sum(DownSpiralSino(:,:,SumSeveralFrames),3);
    sum_mask = sum(DownMask(:,:,SumSeveralFrames),3);
    
    ref_sino = sum_down_cart./sum_mask;
    ref_sino(isnan(ref_sino)) = 0;
    ref_sino(isinf(ref_sino)) = 0;
    
    % % change inverse function by your measurements
    ref = fftshift(fftshift(ifft2(ifftshift(ifftshift(ref_sino,1),2)),2),1);
    
    figure(92);
    imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
    title('Reference image for ME/MC using several frames');
    pause(1);
%     close(92);
    
elseif strcmp(Ref_Type,'FullSingleFrame') == 1
    
    Ref_frame = 1; % define
    
    % % full sampling
    [ref_sino full_mask] = spiral2cart(orig(:,:,Ref_frame),pos(:,:,Ref_frame),nx,ny);
    
    % % change inverse function by your measurements
    ref = fftshift(fftshift(ifft2(ifftshift(ifftshift(ref_sino,1),2)),2),1);
    
    figure(92);
    imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
    title('Reference image for ME/MC using a full single frame');
    pause(1);
    close(92);
end

%% save
save Y.mat DownSpiralSino DownMask Data Pos ref
