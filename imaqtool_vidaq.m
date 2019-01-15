%% Imaqtool Use without Toolbox Attempt
% Plug in Arduino and Hauppage before opening MATLAB

% acquire Device Info for line 9
% info = imaqhwinfo('winvideo');
% ID = info.DeviceInfo(1).DefaultFormat;

close all; clear all; clc;

% Create video input object. 
vid = videoinput('winvideo',1,'UYVY_720x480');

triggerconfig(vid, 'Manual');
vid.FramesPerTrigger = 1;

start(vid)

vid.TriggerRepeat = 9;

tic
for i = 1:10
frame{i} = trigger(vid);
toc
end


stop(vid)