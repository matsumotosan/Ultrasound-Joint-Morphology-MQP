% Script to complete 3D reconstruction (pixel nearest neighbor algorithm)
% of images obtained from 2D ultrasound scans with pose data obtained by an 
% inertial measurement unit. The overall process is as follows:
%
% 0) Enter names of necessary files
% 1) Calculate probe pose from IMU data (filtering, numerical integration)
% 2) US video pre-processing (cropping, pose tagging)
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

% Filter pose data - low pass filter and moving window
% (Filtering function) - ROSIE


% Calculate displacement from linear acceleration using composite
% trapezoidal numerical integration scheme
[disp,vel] = calcDisp(acc,t);


%% 2) US video pre-processing
vid_file = '';  % Enter name of US video file
vid = VideoReader(vid_file);

% Calculate dimension of pixel and cropping window in mm - ROSIE
rect = [];      % cropping window
scale = [];     % mm / pixel

% Convert all frames to grayscale and crop
frames = {};
while hasFrame(vid)
    frames{end + 1} = im2double(imcrop(rgb2gray(readFrame(vid)),rect));
end

% Tag frames with pose data with linear interpolation
frame_pose = interp1(t,pose(:,2:end),frame_time);

%% 3) Calculate bin dimensions and transform coordinates
% Transform IMU pose data to global coordinates - OMEL

% Initialize IMU relatively close to the knee
for i = 1:size(frames,2)
    global_pose = disp .* scale;
end

% Calculate spatial and index limits of bin - OMEL
% I'm thinking about whether this should be global -> local or local ->
% global. I think it would make more sense to go global -> local (mm to px)
% and then convert to real world measurements after.

% Initialize spatial limits of bin
bin_lim = ones(1,6);    % [xmin xmax ymin ymax zmin zmax]

for i = 1:length(pose)
    
    % 3D affine transformation output limits
    tform = affine3d();
    [xlim, ylim, zlim] = outputLimits(tform,[1 bin_lim(2)],[1 bin_lim(4)],[1 bin_lim(6)]);
    
    % Update spatial information
    if xlim(1) < bin_lim(1)
        bin_lim(1) = xlim(1);
    end
    if xlim(2) > bin_lim(2)
        bin_lim(2) = xlim(2);
    end
    if ylim(1) < bin_lim(3)
        bin_lim(3) = ylim(1);
    end
    if ylim(2) > bin_lim(4)
        bin_lim(4) = ylim(2);
    end
    if zlim(1) < bin_lim(5)
        bin_lim(5) = zlim(1);
    end
    if zlim(2) > bin_lim(6)
        bin_lim(6) = zlim(2);
    end
end

% Initialze bin for distribution step
bin = zeros(bin_lim(2), bin_lim(4), bin_lim(6));


%% 4) Reconstruction (distribution step) - SHION
% Fill voxels - PNN algorithm
for i = 1:length(frames)
    bin = fillbin(frames{i},pose(i,:),bin);
    fprintf('(%d/5) frames completed ...\n',i);
end


%% 5) Reconstruction (hole-filling step) - SHION
n = 3;  % grid size
bin = fillhole(bin,n);


%% 6) Visualization


