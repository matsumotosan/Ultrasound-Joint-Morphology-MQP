function bin = fillbin(frame,pose,bin)
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

%% Read bin and pose data
voxSz = size(bin);  % dimension of bin
vox2add = zeros(voxSz(1),voxSz(2),voxSz(3));    % initialize bin

% Translation and rotation values
rx = pose(2);
ry = pose(3);
rz = pose(4);
Sx = pose(5);
Sy = pose(6);
Sz = pose(7);

% Add frame to bin
[fH,fW] = size(frame);
vox2add(1:fH,1:fW,floor(voxSz(3) / 2)) = vox2add(1:fH,1:fW,floor(voxSz(3) / 2)) + frame;


%% Option 1: imwarp
% Affine combinations for rotation and translation

T = eye(4); % 4 by 4 identity matrix
T(1:3,4) = T(1:3,4) + [Sx; Sy; Sz]; % translation in 3D

% Rotation matrices
Rx = [1 0 0 0;
    0 cos(rx) -sin(rx) 0;
    0 sin(rx) cos(rx) 0;
    0 0 0 1];
Ry = [cos(ry) 0 sin(ry) 0;
    0 1 0 0;
    -sin(ry) 0 cos(ry) 0;
    0 0 0 1];
Rz = [cos(rz) -sin(rz) 0 0;
    sin(rz) cos(rz) 0 0;
    0 0 1 0;
    0 0 0 1];

tform = affine3d(Rx * Ry * Rz); % affine 3D transformation object
bin = imwarp(vox2add,tform,'nearest','FillValues',0);

% Compare original and transformed frame
% idx = find(vox2add);
% [a,b,c] = ind2sub(size(vox2add),idx);
% scatter3(a,b,c,vox2add(idx)); hold on
% 
idx = find(bin);
[a,b,c] = ind2sub(size(bin),idx);
scatter3(a,b,c,bin(idx)); hold on
% 
% title('Frame Transformation')
% xlabel('x'); ylabel('y'); zlabel('z');
% legend('Original','Transformed')

%% Option 2: imrotate3
% Fill voxels in output volume outside of limits of original volume with -1
% vox2add = imrotate3(vox2add, pose(2), [1 0 0], 'nearest', 'loose', ...
%     'FillValues', -1);
% vox2add = imrotate3(vox2add, pose(3), [0 1 0], 'nearest', 'loose', ...
%     'FillValues', -1);
% vox2add = imrotate3(vox2add, pose(4), [0 0 1], 'nearest', 'loose', ...
%     'FillValues', -1);


%% Find indices of -1 in output volume
% idx = find(vox2add == -1);
% [I,J,K] = ind2sub(size(vox2add),idx);
% I = unique(I);
% J = unique(J);
% K = unique(K);
% 
% % Update size of vox as necessary
% if (a > voxSz(1))
% %     padarray(vox,[],);
% end
% if (b > voxSz(2))
%     
% end
% if (c > voxSz(3))
%     
% end
% 
% newA = [A, zeros(size(A, 1), size(B, 2)-size(A, 2)); zeros(size(B, 1)-size(A, 1), size(B, 2))];
% vox = vox + vox2add;

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

