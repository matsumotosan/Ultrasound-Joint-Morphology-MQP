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


%% Calculate all index permutations for neighbors
numNeigh = n ^ 3 - 1;   % number of neighbors
shift = -floor(n / 2):floor(n / 2); % index shifts

nbi = zeros(numNeigh,n);    % intialize index shift matrix
di = n;

while di 
    N = n ^ (di - 1);
    ni = 1:N;
    
    while (ni(end) < numNeigh + 1)
        for val = shift
            nbi(ni,di) = val;
            ni = ni+ N;
        end
    end
    
    di = di-1;
end


%% Add the input index array to nbi to get all neighbours
% Set up array for neighbour indices
nd = length(index);
Iadj = zeros(nd,numNeigh);
% kdo = notorig(ones(nd,1), : ); 

for di = 1:ndimA
    indices = imat( :, di );
    shifts = nbi( :, di )';
    neighbindices = indices( :, ones( 1,nneighb)) +shifts( ones(nd, 1), : ) ;
    maxmat = sizeA( di );
    
    % set up mask matrix to keep indices within limits and exclude the original point
    s = logical( neighbindices <= maxmat );
    s =logical( neighbindices > 0 ) & s;
    kdo = kdo & s;
    % Calculate the linear index
    if di == 1       
        Iadj( kdo ) =  neighbindices( kdo );
    else
        blocksize = prod( sizeA( 1:(di-1)  ) );
        m = neighbindices-1;
        Iadj(kdo )  = Iadj(kdo )+ m(kdo)*blocksize;
    end
end

end

