% Script to test pixel nearest neighbor (PNN) hole-filling algorithm. The 
% algorithm scans a pre-filled bin for empty voxels. For each empty voxel,
% the voxel values of the surrounding n-by-n-by-n grid are averaged and
% assigned to the empty voxel.
%
% This test initializes an n-by-n-by-n matrix of zeroes and fills the 
% center of the matrix with ones. The output of the algorithm is 
% demonstrated by comparing 3D scatter plots of the matrix before and 
% after hole-filling.

clear; close all; clc

%% Test hole-filling algorithm (before and after)
% Initialize input matrix
out = 10;  % dimension of zeros matrix
in = 4;    % dimension of inner ones matrix

bin = zeros(out,out,out);   % initialize zeros matrix
idx = (out - in) / 2 + 1:(out + in) / 2;    % ones index vector
[X,Y,Z] = meshgrid(idx,idx,idx);        % indices of ones
% bin(X,Y,Z) = 3 * rand(length(idx) ^ 3 * ones(1,3)); % create inner ones matrix
bin(X,Y,Z) = 100 * ones(length(X)^3,length(X)^3,length(X)^3); % create inner ones matrix

% Hole fill
n = 3;
bin_holefill = fillhole(bin,n);

% Plot before hole-filling
% scatter3(X(:),Y(:),Z(:),40,100 * ones(length(X)^3,1),'filled'); hold on
p = patch(isosurface(1:out,1:out,1:out,bin,100));
isonormals(1:out,1:out,1:out,bin,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 

% Plot after hole-filling
filled_lidx = find(bin_holefill);
[A,B,C] = ind2sub(size(bin_holefill),filled_lidx);
% scatter3(A,B,C,40,bin_holefill(filled_lidx),'filled');



title('Hole-filling Demonstration')
xlabel('x'); ylabel('y'); zlabel('z');
xlim([min(A) max(A)])
ylim([min(B) max(B)])
zlim([min(C) max(C)])
xticks(min(A):max(A))
yticks(min(B):max(B))
zticks(min(C):max(C))

colormap(jet)
colorbar


