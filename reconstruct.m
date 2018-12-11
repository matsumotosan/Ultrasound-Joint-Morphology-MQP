clear; clc; close all

%% IMU data
% IMU reading and settings
run1_imu = cell2mat(struct2cell(load('recon_run1.mat')));
run2_imu = cell2mat(struct2cell(load('recon_run2.mat')));

rate = 0; % sampling rate

% Plot IMU data
t = 0:rate:rate * length(run1_imu);
t = linspace(0,100,length(run1_imu));
plot(t,run1_imu(:,1),t,run1_imu(:,2),t,run1_imu(:,3))

%% Video data
% Load video
dir = '/Users/Shion/Dropbox/MQP-US Media/12-11 test w imu';  % Shion's directory
vid_name = 'recon_run2.mp4';

v = VideoReader(strcat(dir,'/',vid_name));
% noFrames = v.NumberOfFrames;    % number of frames
vHeight = v.Height;
vWidth = v.Width;
vDuration = v.Duration;

% Read frames
frames = {};
i = 1;
rect = [30, 87, 510, 250];

while hasFrame(v)
    frames{i} = imcrop(rgb2gray(readFrame(v)),rect);
%     imshow(frame, 'Parent', currAxes);
%     currAxes.Visible = 'off';
%     pause(0.01);
    i = i + 1;
end


%% Synchrnoize IMU and video



%% Orient frame in 3D
% Initialize 3D matrix of zeros;
vol = zeros(vHeight/2,vHeight,vHeight/2);

% Calculate window height and width from cropping window
crHeight = rect(2) - rect(1);
crWidth = rect(3) - rect(4);


% Orient individual frames in 3D
for i = 1:length(frames)
    f = zeros(vHeight/2,vHeight);
    f = f(
    
    r = zeros(vHeight/2,vHeight,vHeight/2);
    
    r = r(:, :, vHeight/4) + frames{i});
end



