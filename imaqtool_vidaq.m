%% Acquire US video through Image Acquisition Toolbox programmatically
% Plug in Arduino and Hauppage before opening MATLAB
%
% Requires Image Acquisition Toolbox

close all; clear; clc;

%% Acquire images
% acquire Device Info for line 9
% info = imaqhwinfo('winvideo');
% ID = info.DeviceInfo(1).DefaultFormat;

vid = videoinput('winvideo',1,'UYVY_720x480');

triggerconfig(vid, 'Manual');
vid.FramesPerTrigger = 1;

start(vid)
vid.TriggerRepeat = 9;

tic     % start timer
for i = 1:10
    frame{i} = trigger(vid);    % acquire frame
    toc     % time of each acquisition
end

stop(vid)

