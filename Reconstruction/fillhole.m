function bin = fillhole(bin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[x,y,z] = size(bin);    % volume size

% Repeat for each voxel
for i = 1:x
    for j = 1:y
        for k = 1:z
            value = 0;
            count = 0;
            
            if bin(i,j,k) == 0  % voxel is empty
                value = value + 
                
            end
        end
        
        % Assign new voxel value
        bin(i,j,k) = value / counter;
        
    end
end

end

