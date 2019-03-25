% Script to complete 3D reconstruction (pixel nearest neighbor algorithm)
% of images obtained from 2D ultrasound scans with pose data obtained by an 
% inertial measurement unit. The overall process is as follows:
%
% 0) Enter names of necessary files
% 1) Calculate probe pose from IMU data (filtering, numerical integration)
% 2) US video pre-processing (cropping, pose tagging)
% 3) Convert Euler angle to frame rotation
% 3) Calculate bin dimensions and transform coordinates
% 4) Reconstruction (distribution step)
% 5) Reconstruction (hole-filling step)
% 6) Visualization
%
% Requires following MATLAB Toolboxes:
%   - Image Processing
%
% Check toolboxes needed for this script if needed:
% [fList,pList] = matlab.codetools.requiredFilesAndProducts('reconstruct.m');

clear; clc; close all

% Note:
% Try to vectorize as much of the code as possible
% The script is going to be handling a lot of data, so either use MATLAB's
% vectorization or create functions to minimize the number of variables
% needed in the workspace.


%% 0) Enter names of necessary files - comment during development
% imu_file = '';  % raw IMU data
% vid_file = '';  % raw US data

%% 1) Calculate probe pose from IMU data
% Load IMU data
imu_file = 'typr.mat';  % Enter name of IMU data file
imu_data = cell2mat(struct2cell(load(which(imu_file))));

% Filter orientation data
% (Filtering function) - ROSIE


%% 2) US video pre-processing
vid_file = 'recon_run1.mp4';  % Enter name of US video file
addpath(genpath('/Users/Shion/Dropbox/MQP-US Media/3_20 rail tests'))
vid = VideoReader(vid_file);

% Calculate dimension of pixel and cropping window in mm - ROSIE
rect = [29.5 86.5 513 250];      % cropping window
scale = [];     % mm / pixel

% Convert all frames to grayscale and crop
frames = {};

% For video files (.avi, .mp4, etc.)
while hasFrame(vid)
%     frames{end + 1} = im2double(imcrop(rgb2gray(readFrame(vid)),rect));
    frames{end + 1} = double(imcrop(rgb2gray(readFrame(vid)),rect));
end

% Tag frames with pose data with linear interpolation
% pose = interp1(t,pose(:,2:end),frame_time);


%% 3) Calculate bin dimensions

% Calculate mm/pixel
% depth (cm) - cm/pixel (v/h)
% 1 - 0.008/0.010
% 2 - 0.008/0.010
% 3 - 0.008/0.010
% 4 - 0.011/0.014
% 5 - 0.014/0.017
% 6 - 0.017/0.021
% 7 - 0.019/0.024
% 8 - 0.022/0.028
% 9 - 0.025/0.031
% 10 - 
% 11 - 
% 12 - 
% 13 - 
% 14 - 
% 15 - 
% 16 - 

depth = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
mm_per_pixel = {};
res = containers.Map(freq,mm_per_pixel);

%% 4a) Reconstruction (distribution step) - pitch
% Fill voxels - PNN algorithm
% for i = 1:9
%     frames{i} = zeros(15,30);
%     frames{i}(5:10,10:20) = frames{i}(5:10,10:20) + 1;
% end

frames = frames(1:3);
pose = linspace(-45,45,length(frames));

[fH,fW] = size(frames{1});
bin = zeros(2*fH,fW,2*fH,'uint8');

% Fill bin
frames{1} = ones(50,30,'uint8');
for i = 1:8
    frames{end + 1} = 32 * i * ones(50,30,'uint8');
end

bin_ds = fillbin_thick(frames,pose,50,1,1);

% figure
% subplot(1,2,1)
idx = find(bin_ds);
[a,b,c] = ind2sub(size(bin_ds),idx);
% v = 0.7;
% isosurface(a,b,c,bin_ds(idx),v);
scatter3(a,b,c,20,bin_ds(idx),'filled');
colormap(jet)
colorbar
% title('DS')

%% Reconstruction - yaw
noFrames = 3;
angle = linspace(-15,15,noFrames);

bin_thin = fillbin_thin(frames(1:noFrames),angle,300,'nearest');



%% 5a) Reconstruction (hole-filling step)
% n = 9;  % grid size
% bin_hf = fillhole(bin_ds,n);
% 
% subplot(1,2,2)
% idx = find(bin_hf);
% [a,b,c] = ind2sub(size(bin_hf),idx);
% scatter3(a,b,c,20,bin_hf(idx),'filled');
% colormap(jet)
% colorbar
% title('HF')

%% 5b) Reconstruction (3d interpolation)
% tic
% binSz = size(bin_ds);
% idx = find(bin_ds);
% [x,y,z] = ind2sub(size(bin_ds),idx);
% X = [x,y,z];
% v = bin_ds(idx);
% [xx,yy,zz] = meshgrid(1:10:binSz(1),1:10:binSz(2),1:10:binSz(3));
% xq = [xx(:),yy(:),zz(:)];
% vq = griddatan(X,v,xq,'nearest');
% 
% vq = reshape(vq,size(xx));
% slice(xx,yy,zz,vq,80:40:240,100, 100);
% toc

%% 6) Visualization

