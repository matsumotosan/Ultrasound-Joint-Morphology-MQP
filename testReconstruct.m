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
    frames{i} = imcrop(rgb2gray(readFrame(v)),rect);
    i = i + 1;
end

%% Load pose data
pose = [10000 30 45 60 5 10 5];

%% Initialize voxels
vox = zeros(100,100,100,'uint8');
vox_x = linspace(0,300,100);
vox_y = linspace(0,520,100);
vox_z = linspace(0,300,100);
[X,Y,Z] = meshgrid(vox_x,vox_y,vox_z);
voxCoord = [X(:),Y(:),Z(:)];

%% Distribution step
filledBin = fillbin(frames{1}, pose, vox, voxCoord);

%% Hole-filling step


