clear; close all; clc

%% Load video data
dir = '/Users/Shion/Dropbox/MQP-US Media/12-11 test w imu';  % Shion's directory
vid_name = 'recon_run2.mp4';

v = VideoReader(strcat(dir,'/',vid_name));
% noFrames = v.NumberOfFrames;    % number of frames
vDuration = v.Duration;

% Read frames
frames = {};
i = 1;
rect = [30, 87, 510, 250];

while hasFrame(v)
    frames{i} = im2double(imcrop(rgb2gray(readFrame(v)),rect));
    i = i + 1;
end

%% Load pose data
pose = [10000 10 0 0 0 0 0]; % [yaw,pitch,roll,x,y,z]

%% Initialize voxels
% Initialize voxels to hold reconstruction data
vox_dim = [251 511 50];   % reconstruction volume size
vox = zeros(vox_dim(1),vox_dim(2),vox_dim(3));   % initialize matrix of zeros

% Real space 3D coordinates of each voxel
vox_x = linspace(-50,50,vox_dim(1));     % all x-coordinates
vox_y = linspace(-100,100,vox_dim(2));     % all y-coordinates
vox_z = linspace(-50,50,vox_dim(3));     % all z-coordinates

% Create meshgrid to store all coordinate data
[X,Y,Z] = meshgrid(vox_x,vox_y,vox_z);
voxCoord = [X(:),Y(:),Z(:)];

%% Distribution step
filledBin = fillbin(frames{1}, pose, vox, voxCoord);

idx = find(filledBin == 1);
[a,b,c] = ind2sub(vox_dim,idx);
scatter3(a,b,c,'k.'); hold on
xlabel('x')
ylabel('y')
zlabel('z')

%% Hole-filling step


