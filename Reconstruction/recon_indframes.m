clear; clc; close all

%% Read files and create video
folder = '/Users/Shion/Dropbox/MQP-US Media/3_20 rail tests';
addpath(genpath(folder));
file_pattern = fullfile(folder, 'frame_*.mat');
file_list = dir(file_pattern);
frames = {};

% Define cropping window
load('us_setting_map.mat')
depth = 40; % depth setting (mm)
setting = M(depth);
mm_per_pixel = setting(5);
rect = setting(1:4);   % [xmin xmax ymin ymax]
rect = [rect(1), rect(3), rect(2) - rect(1), rect(4) - rect(3)];    % [xmin ymin width height]

% Read and store every image (crop, rgb2gray)
for i = 1:length(file_list)
    imfile = load(file_list(i).name);
    im = imfile.frames;
    frames{i} = double(imcrop(rgb2gray(im{1}), rect));
end


%% Reconstruct (yaw)
% Define probe position
theta_0 = 0;
theta_f = 35;
angles = linspace(theta_0, theta_f, length(frames));
r = 30:40;

for i = 1:length(r)
    % Distribution step
    bin_yaw = fillbin_yaw(frames,angles,r(i),mm_per_pixel,'nearest');

    % Plot
    figure;
    im_min = min(min(bin_yaw(bin_yaw > 0)));    % minimum pixel value
    im_max = max(max(bin_yaw(bin_yaw > 0)));    % maximum pixel value
    bin_yaw(bin_yaw > 0) = bin_yaw(bin_yaw > 0) - im_min;
    imagesc(bin_yaw)
    title('Reconstructed Image')
    xlabel('X (mm)')
    ylabel('Y (mm)')
end




