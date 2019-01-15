function k = findVox(pixCoord,voxCoord,varargin)
%FINDVOX Find voxel index that correspond to coordinates of pixel
%   Return index or indices of voxel with center coordinates closest to
%   provided pixel coordinate. Default is to return single voxel index.
%   Function can also return indices of n-by-n-by-n grid surrounding voxel
%   closest to pixel coordinate.


switch nargin
    case 2
        % N-D nearest point search
        k = dsearchn(voxCoord,pixCoord);
    case 3
        % Dimension of grid surrounding central voxel to return
        n = varargin;
        if ((mod(n,2) == 0) || (floor(n) ~= n) || (n < 0))
            error('n must be odd integer greater than 1')
        end
        xtra = (n - 1) / 2; % number of extra voxels on each side
        a = k(1) - xtra:1:k(1) + xtra;
        b = k(2) - xtra:1:k(2) + xtra;
        c = k(3) - xtra:1:k(3) + xtra;

        % Indices of surrounding n-by-n-by-n grid
        [A,B,C] = meshgrid(a,b,c);
        k = [A(:)'; B(:)'; C(:)'];
    otherwise
        error('Function only accepts 2 or 3 arguments')
end

end

