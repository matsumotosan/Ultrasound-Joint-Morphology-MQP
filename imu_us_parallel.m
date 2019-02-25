%% Acquire IMU data and US images simultaneously
% Requires Parallel Processing Toolbox
% This script contains program to collect US images

% Close vid instance if one exists
if exist('vid')
    stop(vid);
end

clear; close all; clc

%% Idea 1: Open another MATLAB instance
% Evaluate terminal command through MATLAB
mr = matlabroot;            % location of MATLAB in directory
sr = which('gyroviz.m');    % path to gyroviz.m

% Windows - evaluate command in bash
eval(strcat('!matlab -nodesktop -nosplash -r'," ", '"run(', "'", sr, "'",')"'))
% eval('!matlab -nodesktop -nosplash -r "gyroviz.m" &')

% Allow communication between scripts (need TCP/UDP/IP Toolbox 2.0.6)


% Wait for IMU to stabilize


%% Idea 2: Parallel processing
% parpool
% 
% spmd
%     
% end

if exist('vid')
    stop(vid);
end

% Acquire images
vid = videoinput('winvideo',1,'UYVY_720x480');
preview(vid);

triggerconfig(vid, 'Manual');

delete(gcp('nocreate'))
parpool(2)
parfor idx = 1:2
    if idx == 1
        test1();          %IMU data acquisition
    elseif idx == 2
        test2();   %image acquisition from ultrasound
    end
end

%% Acquire images
% Acquire Device Info for line 9
% info = imaqhwinfo('winvideo');
% ID = info.DeviceInfo(1).DefaultFormat;

% Create video object
vid = videoinput('winvideo',1,'UYVY_720x480');
preview(vid);   % Preview video

% Configure video object to only acquire frame when triggered
triggerconfig(vid, 'Manual');

% Gather Frames
start(vid)
tic

noFrames = 4000;    % number of frames to be acquired
frames = {};        % initialize cell to hold frames

for i = 1:noFrames
    
    frames{i} = getsnapshot(vid);
    t(i) = toc;

    % potential timestamp tag
    %     if i == 1
    %         fprintf('Start time of data acquistion: %s \n', datestr(now,'mm-dd-yyyy MM.SS.FFF'));
    %     end
end

% Close Video
stop(vid)
delete(vid)


