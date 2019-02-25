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


%% Initialize voxels - frames from US
% Load pose data
pose = [10000 10 0 0 0 0 0]; % [time,yaw,pitch,roll,Sx,Sy,Sz]
Sx = pose(:,5);
Sy = pose(:,6);
Sz = pose(:,7);

% Frame data
noFrames = length(frames);
[fH,fW] = size(frames{1});
c = 0.1;    % mm / pixel

% Calculate dimensions of voxel to initialize
H = ceil((max(Sy) - min(Sy)) / c) + fH;  % height
W = ceil((max(Sx) - min(Sx)) / c) + fW;  % width
D = ceil((max(Sz) - min(Sz)) / c) + fW;  % depth

% Calculate output limits for 3D affline transformation
% [xlim,ylim,zlim] = outputLimits(tform,x,y,z);

% Initialize voxels to hold reconstruction data
% vox = zeros(H,W,D);   % initialize matrix of zeros
vox = zeros(300, 600, 100);

% Real space 3D coordinates of each voxel
vox_x = linspace(-50,50,W);     % all x-coordinates
vox_y = linspace(-100,100,H);   % all y-coordinates
vox_z = linspace(-50,50,D);     % all z-coordinates

% Create meshgrid to store all coordinate data
[X,Y,Z] = meshgrid(vox_x,vox_y,vox_z);
voxCoord = [X(:),Y(:),Z(:)];

%% Verification cube
cube = zeros(100,100,100);
cube(25:75,25:75,25:75) = cube(25:75,25:75,25:75) + 1;


%% Distribution step
clear frames
for i = 1:5
    frames{i} = rand(20,20);
end

filledBin = zeros(300,300,300);
pose = [10000 0 0 0 0 0 0;
    10000 0 0 0 0 0 0;
    10000 0 0 0 0 0 0;
    10000 0 0 0 0 0 0;
    10000 0 0 0 0 0 0];

%tic
for i = 1:4
    filledBin = fillbin(frames{i}, pose(i,:), filledBin);
    fprintf('(%d/5) frames completed ...\n',i);
end
%toc

% idx = find(filledBin);
% [a,b,c] = ind2sub(size(filledBin),idx);
% scatter3(a,b,c,filledBin(idx));
% 
% title('Reconstruction')
% xlabel('x'); ylabel('y'); zlabel('z');

%% Hole-filling step


