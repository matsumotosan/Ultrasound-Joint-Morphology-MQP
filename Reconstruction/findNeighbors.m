function neighbors = findNeighbors(index,size,n)
%FINDNEIGHBORS Find linear indices of voxels surrounding index
%   Algorithm adopted from neighborND function written by Ronald Ouwekerk.
%   Extended to find linear indices of neighbor points in an n-by-n-by-n
%   grid. Remove spatial resolution functionality.
%
% Input:  index = subscript index of center voxel
%          size = bin dimensions
%             n = grid size surrounding center voxel  
%
% Output:   idx = linear indices of neighboring voxels

% n must be an odd integer greater than 1
if ~mod(n,2) || (n <= 1)
    error('n must be an odd integer greater than one');
end


%% Calculate all index permutations
shift = -floor(n / 2):floor(n / 2);
nbi = permn(shift,3);
nbi(~any(nbi,2),:) = [];


%% Calculate indices of neighbors
% Add shift permutations to center voxel index
neighbors = index + nbi;

% Remove invalid indices (outside bounding box)
for i = 1:3
    neighbors(any(neighbors(:,i) > size(i),2),:) = [];
    neighbors(any(neighbors(:,i) < 1,2),:) = [];
end


%% Plot results (comment for speed)
% Scatter plot of neighbors and center voxel
% scatter3(neighbors(:,1),neighbors(:,2),neighbors(:,3),'MarkerFaceColor',[0.75 0 .75]); hold on
% scatter3(index(1),index(2),index(3),'MarkerFaceColor',[0 .75 .75]);
% 
% % Plot representative voxels around all points
% cb = 0.3 * ones(1,3);
% 
% % figure; hold on
% plotcube(cb,index - 0.5 * cb(1),0.5,[1 0 0]); hold on % center voxel
% for i = 1:length(neighbors) % neighbor voxels
%     plotcube(cb,neighbors(i,:) - 0.5 * cb(1),0.5,[0 1 0]);
% end
% 
% grid on
% % legend('Center Voxel','Neighboring Voxels');
% title('Voxel Neighbor Finding Algorithm');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% xlim([min(neighbors(:,1)) - cb(1) * 0.5, max(neighbors(:,1)) + cb(1) * 0.5])
% xlim([min(neighbors(:,2)) - cb(1) * 0.5, max(neighbors(:,2)) + cb(1) * 0.5])
% xlim([min(neighbors(:,3)) - cb(1) * 0.5, max(neighbors(:,3)) + cb(1) * 0.5])
% xticks(min(neighbors(:,1)):max(neighbors(:,1)))
% yticks(min(neighbors(:,2)):max(neighbors(:,2)))
% zticks(min(neighbors(:,3)):max(neighbors(:,3)))

end

