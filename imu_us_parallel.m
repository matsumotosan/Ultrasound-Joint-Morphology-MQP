%% Acquire IMU data and US images simultaneously
% Requires Parallel Processing Toolbox
% This script contains program to collect US images
clear; close all;

% Close vid instance if one exists
if exist('vid')
    stop(vid);
end

clc

%% Idea 1: Open another MATLAB instance
% Evaluate terminal command through MATLAB
mr = matlabroot;    % location of MATLAB in directo

% Windows
eval('!matlab -nodesktop -nosplash -r "gyroviz.m" &')


% MAC
% Add following to end of .bash_profile
% (located in home directory - type "cd ~/.bash_profile" in Terminal)
%
% export PATH=$PATH:/Applications/MATLAB_R2015b.app/bin/
% eval('!osascript -e ''tell application "Terminal"'' -e ''activate'' -e ''do script "matlab -nodesktop -nosplash -r \"gyroviz.m\""'' -e ''end tell''')

% Allow communication between scripts (need TCP/UDP/IP Toolbox 2.0.6)


% Wait for IMU to stabilize


%% Idea 2: Parallel processing
% parpool
% 
% spmd
%     
% end

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


