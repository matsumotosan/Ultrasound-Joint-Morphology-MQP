function us_data = collect_us(noframes)
%% Gather Frames
start(vid)

us_data = {};

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