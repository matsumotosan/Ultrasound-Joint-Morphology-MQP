function bin = fillbin(frame,pose,bin,r,opt)
%FILLBIN Traverse input pixels of frame and insert them into corresponding
%voxels
%   Distribution (DS) step of pixel-based method known as pixel nearest
%   neighbor (PNN) method. If multiple pixels correspond to one voxel, the
%   pixel values are averaged.
%
% Input:     frame = grayscale 2D scan
%             pose = pose of frame (tilt)
%              vox = matrix containing all voxel data
%         voxCoord = coordinate of each voxel
%
% Output:      vox = matrix of voxels with frame added according to PNN
%                    method
%

if strcmp(opt,'yaw')
    axis = [1 0 0];
elseif strcmp(opt,'pitch')
    axis = [0 1 0];
else
    error('opt must be either yaw or pitch')
end

%% Read bin and pose data
% Frame and bin data
voxSz = size(bin);              % bin dimensions
[fH,~] = size(frame{1});        % frame height
defX = round(voxSz(1) / 2);     % default frame x-coordinate
defZ = voxSz(3) - fH + 1:voxSz(3);  % default frame z-coordinates

% Add each frame to bin
for i = 1:length(frame)
    % Initialize bin of zeros
    vox2add = zeros(voxSz(1),voxSz(2),voxSz(3));

    % Add frame
    vox2add(defX,:,defZ) = squeeze(vox2add(defX,:,defZ)) + frame{i}';
    
    % Rotate frame
    vox2add = imrotate3(vox2add,pose(i),axis,'nearest','crop');
    
    % Translate frame
    dx = r * sind(abs(pose(i)));
    dz = r * (1 - cosd(abs(pose(i))));
    if pose(i) < 0
        vox2add = imtranslate(vox2add, [0 -dx -dz]);
    elseif pose(i) > 0
        vox2add = imtranslate(vox2add, [0 dx -dz]);
    end
    
    % Add to bin
    bin = bin + vox2add;
end

% Average all voxels
% bin = bin / length(frame);

% Plot result
% idx = find(bin);
% [a,b,c] = ind2sub(size(bin),idx);
% figure
% scatter3(a,b,c); hold on
% 
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% 
% xlim([0 voxSz(1)])
% ylim([0 voxSz(2)])
% zlim([0 voxSz(3)])


end

