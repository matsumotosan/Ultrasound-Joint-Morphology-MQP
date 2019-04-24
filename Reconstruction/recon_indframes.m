clear; clc; close all

%% Read files and create video
% folder = '/Users/Shion/Dropbox/MQP-US Media/4_2_final_test/with_bag/trial_2';
% folder = '/Users/Shion/Dropbox/MQP-US Media/4_6_phantom_test2/trial3/angles_fixed';
folder = '/Users/Shion/Dropbox/MQP-US Media/3_26_rail_test';
addpath(genpath(folder));
file_pattern = fullfile(folder, '*.mat');
file_list = dir(file_pattern);
frames = {};

% Define cropping window
load('us_setting_map.mat')
depth = 50; % depth setting (mm)
setting = M(depth);
mm_per_pixel = setting(5);
rect = setting(1:4);   % [xmin xmax ymin ymax]
rect = [rect(1), rect(3), rect(2) - rect(1), rect(4) - rect(3)];    % [xmin ymin width height]

% Read and store every image (crop, rgb2gray) and its angle
theta = zeros(1,length(file_list));
for i = 1:length(file_list)
    imfile = load(file_list(i).name);
    im = imfile.frames;
    figure
    imagesc(im{1})
    frames{i} = double(imcrop(rgb2gray(im{1}), rect));
    [~,name,~] = fileparts(file_list(i).name);
    C = strsplit(name,'_');
    theta(i) = str2double(C{4});
end


%% Reconstruct (yaw)
% Define probe position
theta = theta - min(theta);
r_nom = 40;
r = linspace(r_nom - 5, r_nom + 5, 12);  % rail imaging radius (mm)
r = r_nom;
method = {'nearest','bilinear','bicubic'};

% figure;
for i = 1:length(r)
    % Distribution step
    bin_yaw = fillbin_yaw(frames(:),theta,r(i),mm_per_pixel,method{1});

    % Change scaling so minimum pixel value is zero
    im_min = min(min(bin_yaw(bin_yaw > 0)));    % minimum pixel value greater than 0
    im_max = max(max(bin_yaw));                 % maximum pixel value
    bin_yaw(bin_yaw > 0) = bin_yaw(bin_yaw > 0) - im_min;
    
    % Plot reconstructed image
    figure
    imagesc(bin_yaw)
    title(['Reconstructed Image (r=' num2str(r(i)) 'mm)'])
    xlabel('X (mm)')
    ylabel('Y (mm)')
    set(gca,'TickDir','out');
    noTicks = 10;
    
    xticks(linspace(0,size(bin_yaw,2), noTicks));
    xticklabels(round(linspace(0,size(bin_yaw,2) * mm_per_pixel,noTicks),1));
    xtickangle(45)
    yticks(linspace(0,size(bin_yaw,1),noTicks));
    yticklabels(round(linspace(0,size(bin_yaw,1) * mm_per_pixel,noTicks),1));
    
    colorbar
    
    % Create distance tool
    % To calculate distance between two points, multiply the distance in
    % pixel units by mm_per_pixel to get the Euclidean distance in mm
    h = imdistline(gca, [10 100], [10 100]);
end




