function bin = fillbin_thin(frame,angle,r,method)
%FILLBIN_THIN - Distribute pixels to representative voxels

%%

[dims,rows,cols,p] = init_bin(frame{1},r);
bin = zeros(dims(1),dims(2),'double');

for i = 1:length(frame)
    % Initialize bin
    bin2add = zeros(dims(1),dims(2),'double');
    
    % Add frame
    bin2add(rows,cols(1:end-1)) = bin2add(rows,cols(1:end-1)) + frame{i};
    
    % Rotate frame around point
    bin2add = rotateAround(bin2add,p(1),p(2),angle(i),method);
    
    % Add to bin
    bin = bin + bin2add;
end

bin = floor(bin / length(frame));

% Grayscale image
imshow(uint8(bin)); hold on
plot(p(1),p(2),'r+');   % rotation point
title('Grayscale Reconstruction')
xlabel('Horizontal')
ylabel('Depth')

end

function [dims,rows,cols,point] = init_bin(frame,r)

    % Pixel dimension of frames
    [fh,fw] = size(frame);

    % Reconstruction window size
    dims = [fh + r + ceil(fw / 2),2 * (fh + r)];

    % Index of initial image insertion in placeholder bin
    rows = 1:fh;
    cols = floor((dims(2) - fw) / 2:(dims(2) - fw) / 2 + fw); 
    
    % Point to rotate around
    point = [fh + r + 1,dims(2) / 2];
end
