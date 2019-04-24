%% Acquire IMU data and US images simultaneously
% Requires Parallel Processing Toolbox
% This script is intended to collect US scans and IMU data simultaneously
% to align US scans and IMU data temporally.

clear; close all; clc

%% Idea 1: Open another MATLAB instance
% % Evaluate terminal command through MATLAB
% mr = matlabroot;            % location of MATLAB in directory
% sr = which('gyroviz.m');    % path to gyroviz.m
% 
% % Windows - evaluate command in bash
% eval(strcat('!matlab -nodesktop -nosplash -r'," ", '"run(', "'", sr, "'",')"'))
% % eval('!matlab -nodesktop -nosplash -r "gyroviz.m" &')
% 
% % Allow communication between scripts (need TCP/UDP/IP Toolbox 2.0.6)
% 
% 
% % Wait for IMU to stabilize


%% Idea 2: Parallel MATLAB workers
% if exist('vid')
%     stop(vid);
% end

% IMU info
baudrate = 115200;
port = 'COM5';   % Arduino UNO
%port = 'COM6';   % Arduino Nano
n = 1000;        % number of IMU measurements

% US info
noframes = 100;

delete(gcp('nocreate'))
parpool(2)

fname = 'trial1' + '%d.mat';

% Collect IMU and US data in parallel
parfor idx = 1:2
    if idx == 1
        imu_data = collect_imu(baudrate,port,n);    %IMU data acquisition
    elseif idx == 2
        us_data = collect_us(noframes);             % US image acquisiton
    end
    
    % Save variables
    parsave(sprintf(fname,idx),x,y);
end

%% Acquire images
% Acquire Device Info for line 9
% info = imaqhwinfo('winvideo');
% ID = info.DeviceInfo(1).DefaultFormat;

% Create video object
% vid = videoinput('winvideo',1,'UYVY_720x480');
% preview(vid);   % Preview video
% 
% % Configure video object to only acquire frame when triggered
% triggerconfig(vid, 'Manual');
% 
% % Gather Frames
% start(vid)
% tic
% 
% noFrames = 4000;    % number of frames to be acquired
% frames = {};        % initialize cell to hold frames
% 
% for i = 1:noFrames
%     
%     frames{i} = getsnapshot(vid);
%     t(i) = toc;
% 
%     % potential timestamp tag
%     %     if i == 1
%     %         fprintf('Start time of data acquistion: %s \n', datestr(now,'mm-dd-yyyy MM.SS.FFF'));
%     %     end
% end
% 
% % Close Video
% stop(vid)
% delete(vid)

%%
parfor ii = 1:4
    x = rand(10,10);
    y = ones(1,3);
    parsave(sprintf('output%d.mat', ii), x, y);
end
