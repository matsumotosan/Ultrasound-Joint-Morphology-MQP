function bin = fillbin_thin(frame,angle,r,method)
%FILLBIN_THIN - Distribute pixels to representative voxels
%
%   FILLBIN_THIN(FRAME,ANGLE,R,METHOD) - performs distribution step of
%   reconstruction for probe in yaw motion. Each FRAME is inserted in a
%   bin at the corresponding ANGLE on a rail of radius R. Pixel
%   interpolation specified by METHOD:
%       'nearest'       Neares-neighbor interpolation (fastest)
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation (best results)
%
%%

[dims,rows,cols,p] = init_bin(frame{1},r);  % bin info
bin = zeros(dims(1),dims(2),'double');      % initialize bin

for i = 1:length(frame)
    % Initialize bin
    bin2add = zeros(dims(1),dims(2),'double');
    
    % Add frame
    bin2add(rows(1:end-1),cols(1:end-1)) = bin2add(rows(1:end-1),cols(1:end-1)) + frame{i};
    
    % Rotate frame around point
    bin2add = rotateAround(bin2add,p(1),p(2),angle(i),method);
    
    % Add to bin
    bin = bin + bin2add;
end

bin = floor(bin / length(frame));

% Grayscale image
figure; hold on
imshow(uint8(bin)); 
% colormap jet
% colorbar

% plot(p(2),p(1),'r+','MarkerSize',20);   % rotation point
title('Grayscale Reconstruction')
xlabel('Horizontal')
ylabel('Depth')

end

function [dims,rows,cols,point] = init_bin(frame,r)
    % Pixel dimension of frames
    [fh,fw] = size(frame);

    % Reconstruction window size
    dims = [r + fw,2 * (fh + r)];

    % Index of initial image insertion in placeholder bin
    rows = floor(fw / 2:fw / 2 + fh);
    cols = floor((dims(2) - fw) / 2:(dims(2) - fw) / 2 + fw); 
    
    % Point to rotate around
    point = [round(fw + fh) / 2 + r + 1,dims(2) / 2];
end
