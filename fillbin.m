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

voxSz = size(vox);
vox2add = zeros(voxSz(1),voxSz(2),voxSz(3));

vox2add(:,:,floor(voxSz(3) / 2)) = vox2add(:,:,floor(voxSz(3) / 2)) + frame;
vox2add = imrotate3(vox2add, pose(2), [1 0 0]);
vox2add = imrotate3(vox2add, pose(3), [0 1 0]);
vox2add = imrotate3(vox2add, pose(4), [0 0 1]);

% Update size of vox if necessary
[a,b,c] = size(vox2add);
if (a > voxSz(1))
    padarray();
end
if (b > voxSz(2))
    
end
if (c > voxSz(3))
    
end
newA = [A, zeros(size(A, 1), size(B, 2)-size(A, 2)); zeros(size(B, 1)-size(A, 1), size(B, 2))];

vox = vox + vox2add;

% % Coordinates of frame corners before transformation (x1 x2 y1 y2 z)
% defaultCoord = [0 50 0 100 0];
% 
% [h,w] = size(frame);    % size of frame 
% voxSz = size(vox);      % size of vox
% 
% % Find transformed 3D-coordinates of pixels
% pixCoord = findPos(pose,[h w],defaultCoord);
% 
% % Repeat for every pixel in frame
% for i = 1:h
%     for j = 1:w
%         
%         % Find corresponding voxel indices - returns row of voxCoord
%         voxIdx = findVox(pixCoord(w * (i - 1) + j,:), voxCoord);
%         [I,J,K] = ind2sub(voxSz,voxIdx);
%         
%         if vox(I,J,K) == 0  % voxel is empty
%             % Fill voxel with pixel value
%             vox(I,J,K) = frame(i,j);
%             pixCount = 1;
%         else    % voxel has a value
%             % Average with new pixel value
%             vox(I,J,K) = (vox(I,J,K) * pixCount + frame(i,j)) / (pixCount + 1);
%             pixCount = pixCount + 1;    % increment pixel count
%         end
%     end
% end

end

