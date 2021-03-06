function us_data = collect_us(noframes)
% COLLECT_US - acquire NOFRAMES frames of US scans
%
%   COLLECT_US(NOFRAMES) - to be run using imu_us_parallel.m script to
%   acquire US scans in parallel with IMU data acquisition. Acquire
%   NOFRAMES frames before exiting.

%% 
% Initialize cell array
us_data = cell(noframes,2);

% Acquire images
vid = videoinput('winvideo',1,'UYVY_720x480');
preview(vid);

triggerconfig(vid, 'Manual');

start(vid)

tic
for i = 1:noframes
    
    % Gather frame and time data
    us_data{i,1} = getsnapshot(vid);
    us_data{i,2} = toc;

    % potential timestamp tag
    %     if i == 1
    %         fprintf('Start time of data acquistion: %s \n', datestr(now,'mm-dd-yyyy MM.SS.FFF'));
    %     end
end

%% Close Video
stop(vid)
delete(vid)

%% Calculate fs

% t = t';
% fs = 1/mean(diff(t));

end