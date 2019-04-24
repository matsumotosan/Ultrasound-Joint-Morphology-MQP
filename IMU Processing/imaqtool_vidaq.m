%% Acquire US video through Image Acquisition Toolbox
% Plug in Arduino and Hauppage before opening MATLAB
% If you encounter any issues unplug the USB for Hauppauge, close MATLAB,
% plug Hauppauge back in and then open a new MATLAB session.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This code requires Image Acquisition Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For use of the GUI associated with the Image Acq. Toolbox: type the
% following command into the command window:
% imaqtool

%% Reset the MATLAB session
close all; clear; clc;

% delete video object if it exists
if exist('vid')
    stop(vid);
end

%% Acquire images

% % acquire Device Info 1for line 17
% info = imaqhwinfo('winvideo');
% ID = info.DeviceInfo(1).DefaultFormat;

%create video object & preview
vid = videoinput('winvideo',1,'UYVY_720x480'); %inputs may change per computer
preview(vid);

%allow for manual trigger
triggerconfig(vid, 'Manual');

%% Gather Frames
start(vid)

%change frame num
frame_num = input('How many frames do you want? [#]  '); 

%begin timer
tic

%iterate 'getsnapshot' for number of desired frames
for i = 1:frame_num
    
    %save frame
    frames{i} = getsnapshot(vid);
    
    %store time per frame
    t(i) = toc;

    % potential timestamp tag
    %     if i == 1
    %         fprintf('Start time of data acquistion: %s \n', datestr(now,'mm-dd-yyyy MM.SS.FFF'));
    %     end
end

%% Close Video
stop(vid)
delete(vid)

%% Calculate sampling frequency (fs)

%reconfigure time
t = t';

%calculate sampling frequency
fs = 1/mean(diff(t));
