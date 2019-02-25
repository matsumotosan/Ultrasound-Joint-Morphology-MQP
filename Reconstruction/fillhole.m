function bin = fillhole(bin,n)
%FILLHOLE Pixel nearest neighbor hole-filling with averaging algorithm
%   Algorithm adapted from "Freehand 3D Ultrasoud Reconstruction Algorithms -
%   A Review". Authors: Solberg, Lindseth, Torp, et al., 2007
%
% Input:  bin = bin after reconstruction
%           n = n-by-n-by-n search grid for nonzero voxels around empty voxels
%
% Output: bin = bin with PNN hole-filling with averaging
%


% n must be an odd integer greater than 1
if ~mod(n,2) || (n <= 1) || ~isnumeric(n)
    error('n must be an odd integer greater than one');
end

sz = size(bin);              % bin dimensions
lidx = find(~bin);           % indices of empty voxels
[I,J,K] = ind2sub(sz,lidx);  % convert to subscript index

% Repeat for each empty voxel
for i = 1:length(I)
    value = 0;  % sum of neighboring voxel values

    % linear indices of neighboring voxels
    neigh = findNeighbors([I(i),J(i),K(i)],sz,n);
            
    % calculate average of nonzero voxels
    for j = 1:length(neigh)
        if bin(neigh(j))
            value = value + bin(neigh(j));
        end
    end
    
    % Assign average value or zero to voxel
    bin(I(i),J(i),K(i)) = value / j;
end

end

