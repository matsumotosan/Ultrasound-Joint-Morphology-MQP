function vox = fillbin(frame,pose,vox)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[h,w] = size(frame);    % size of frame

% Repeat for every pixel in frame
for i = 1:h
    for j = 1:w
        % Find voxel the pixel belongs to
        [a,b,c] = findVox(pose);
        
        if vox(a,b,c) == 0  % voxel is empty
            % Fill voxel with pixel value
            vox(a,b,c) = frame(i,j);
            pixCount = 1;
        else    % voxel has a value
            % Average with new pixel value
            vox(a,b,c) = (vox(a,b,c) * pixCount + frame(i,j)) / (pixCount + 1);
            pixCount = pixCount + 1;    % increment pixel count
        end
    end
end



end

