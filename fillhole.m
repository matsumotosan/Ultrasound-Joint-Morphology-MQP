function vox = fillhole(vox,grid)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[x,y,z] = size(vox);    % volume size

% Repeat for each voxel
for i = 1:x
    for j = 1:y
        for k = 1:z
            value = 0;
            count = 0;
            
            if vox(i,j,k) == 0  % voxel is empty
                value = value + 
                
            end
        end
        
        % Assign new voxel value
        vox(i,j,k) = value / counter;
        
    end
end

end

