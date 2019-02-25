function test2()
%% Gather Frames
start(vid)

%change frame num
frame_num = 10

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

end