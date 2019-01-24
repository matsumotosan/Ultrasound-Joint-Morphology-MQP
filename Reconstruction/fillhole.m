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
if ~mod(n,2) || (n <= 1)
    error('n must be an odd integer greater than one');
end

sz = size(bin);    % bin dimensions

% Repeat for each voxel
for i = 1:sz(1)
    for j = 1:sz(2)
        for k = 1:sz(3)
            value = 0;
            count = 0;
            
            if ~bin(i,j,k)  % voxel is empty
                neigh = findNeighbors([i,j,k],sz,n);
                for n = length(neigh)
                    
                end
            end
        end
        
        % Assign averaged voxel value
        bin(i,j,k) = value / counter;
    end
end

end

