
%% Add path 
clear all;
addpath('bin');
addpath('data');

%% Downsampling rate 
DsRate = 6;     % define
disp(['Downsample rate: ',num2str(DsRate)]);

%% Select Reference type for ME/MC
Ref_Type = 'FullSingleFrame';
% Ref_Type = 'DiastoleFrames';

%% Read full measurement 
% % Cartesian case (fft(image,[],1))
FileName = ['full_cart.mat']; % file name => sino
load(FileName);

Scaling = 1;    % define
sino = sino*Scaling;

%% you should set diastole frames or a single full frame for making the reference image for ME/MC procedure.

% low frequency full sampling number (1~num_phase/2), (end-num_phase/2~end) -> full sampling
num_low_phase = 8;      % define

[nx ny nt] = size(sino);

M = Random_DownsamplingMASK(nx,ny,nt,DsRate,num_low_phase); % Mask generation

if strcmp(Ref_Type,'DiastoleFrames') == 1
    
    Diastoleframes = [15:25]; % define

    % non-overlap masking for diastole frames
    M(:,:,Diastoleframes) = Random_DownsamplingMASK_unif_nooverlap(nx, ny, length(Diastoleframes), DsRate,num_low_phase);
    
    DownSino = M.*sino;
    mask = M;
    
    sum_down_cart = sum(DownSino(:,:,Diastoleframes),3);
    sum_mask = sum(mask(:,:,Diastoleframes),3);
 
    ref_sino = sum_down_cart./sum_mask;
    ref_sino(isnan(ref_sino)) = 0;
    ref_sino(isinf(ref_sino)) = 0;

    ref = ifft(ref_sino,[],1);

    figure(92);
    imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
    title('Reference image for ME/MC using diastole frames');
    pause(1);
    close(92);
    
elseif strcmp(Ref_Type,'FullSingleFrame') == 1
    
    DownSino = M.*sino;
    mask = M;
    
    Ref_frame = 1; % define 
    
    ref = ifft(sino(:,:,Ref_frame),[],1);

    figure(92);
    imagesc(abs(ref)); axis off; axis equal; colormap gray; colorbar;
    title('Reference image for ME/MC using a full single frame');
    pause(1);
    close(92);
    
end

%% save 
save Y.mat DownSino mask ref num_low_phase
