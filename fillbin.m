function vox = fillbin(frame,pose,vox,voxCoord)
%FILLBIN Traverse input pixels of frame and insert them into corresponding
%voxels
%   Distribution (DS) step of pixel-based method known as pixel nearest
%   neighbor (PNN) method. If multiple pixels correspond to one voxel, the
%   pixel values are averaged.
%
% Input:     frame = grayscale 2D scan
%             pose = pose of frame (yaw, pitch, roll, x-disp, y-disp, z-disp)
%              vox = matrix containing all voxel data
%         voxCoord = coordinate of each voxel
%
% Output:      vox = matrix of voxels with frame added according to PNN
%                    method
%

% Coordinates of frame before transformation (x1 x2 y1 y2 z)
defaultCoord = [0 50 0 100 0];

[h,w] = size(frame);    % size of frame
[xLen,yLen,~] = size(vox);   % size of vox    

% Find transformed 3D-coordinates of pixels
pixCoord = findPos(pose,[h w],defaultCoord);

% Repeat for every pixel in frame
for i = 1:h
    for j = 1:w
        
        % Find corresponding voxel indices - returns row of voxCoord
        voxIdx = findVox(pixCoord(i*j + j - 1,:), voxCoord);
        vox_xIdx = mod(mod(voxIdx / (xLen * yLen)), xLen) + 1;
        vox_yIdx = mod(voxIdx,yLen) + 1;
        vox_zIdx = floor(voxIdx / (xLen * yLen)) + 1;
        
        if vox(vox_xIdx,vox_yIdx,vox_zIdx) == 0  % voxel is empty
            % Fill voxel with pixel value
            vox(vox_xIdx,vox_yIdx,vox_zIdx) = frame(i,j);
            pixCount = 1;
        else    % voxel has a value
            % Average with new pixel value
            vox(vox_xIdx,vox_yIdx,vox_zIdx) = (vox(vox_xIdx,vox_yIdx,vox_zIdx) ...
                * pixCount + frame(i,j)) / (pixCount + 1);
            pixCount = pixCount + 1;    % increment pixel count
        end
    end
end

end

