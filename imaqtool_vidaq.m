%% Acquire US video through Image Acquisition Toolbox programmatically
% Plug in Arduino and Hauppage before opening MATLAB
% Best when running alone on computer
% Requires Image Acquisition Toolbox

close all; clear; clc;

if exist('vid')
    stop(vid);
end

%% Acquire images
% acquire Device Info for line 17
% info = imaqhwinfo('winvideo');
% ID = info.DeviceInfo(1).DefaultFormat;

vid = videoinput('winvideo',1,'UYVY_720x480');
preview(vid);

triggerconfig(vid, 'Manual');

%% Gather Frames
start(vid)

%change frame num
frame_num = input('How many frames do you want? [#]  '); 

tic
for i = 1:frame_num
    
    frames{i} = getsnapshot(vid);
    t(i) = toc;

    % potential timestamp tag
    %     if i == 1
    %         fprintf('Start time of data acquistion: %s \n', datestr(now,'mm-dd-yyyy MM.SS.FFF'));
    %     end
end


%% Close Video
stop(vid)
delete(vid)

%% Calculate fs

t = t';
fs = 1/mean(diff(t));
