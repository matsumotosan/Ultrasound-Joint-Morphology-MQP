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

% Plot original slice
idx = find(vox2add);
[a,b,c] = ind2sub(size(vox2add),idx);
figure
scatter3(a,b,c,10); hold on


%% Option 1: Quaternions



%% Option 2: imwarp (rotation matrix)
% % Rotation matrices
% Rx = [1 0       0        0;
%       0 cos(rx) -sin(rx) 0;
%       0 sin(rx) cos(rx)  0;
%       0 0       0        1];
% Ry = [cos(ry)  0 sin(ry) 0;
%       0        1 0       0;
%       -sin(ry) 0 cos(ry) 0;
%       0        0 0       1];
% Rz = [cos(rz) -sin(rz) 0 0;
%      sin(rz)  cos(rz)  0 0;
%      0        0        1 0;
%      0        0        0 1];
% 
% % Transformation
% tform = affine3d(Rx * Ry * Rz); % affine 3D transformation object
% vox2add = imwarp(vox2add,tform,'nearest','FillValues',0);   % rotate image
% % vox2add = imtranslate(vox2add,[Sx Sy Sz],'nearest','FillValues',0); % translate image
% bin = bin + vox2add;
% 
% % Compare original and transformed frame
% idx = find(vox2add);
% [a,b,c] = ind2sub(size(vox2add),idx);
% scatter3(a,b,c,10,vox2add(idx)); hold on
% 
% idx = find(bin);
% [a,b,c] = ind2sub(size(bin),idx);
% scatter3(a,b,c,10); hold on
% 
% title('Frame Transformation')
% xlabel('x'); ylabel('y'); zlabel('z');
% legend('Original','Transformed')


end

