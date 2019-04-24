function bin = fillbin_yaw(frame,angle,r,mm_per_pixel,method)
%FILLBIN_YAW - Reconstruction in yaw distribution step
%
%   FILLBIN_YAW(FRAME,ANGLE,R,MM_PER_PIXEL,METHOD) - performs distribution 
%   step of reconstruction for probe in yaw motion. Each FRAME is inserted 
%   in a bin at the corresponding ANGLE on a rail of radius R. Pixel
%   interpolation specified by METHOD:
%       'nearest'       Neares-neighbor interpolation (fastest)
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation (best results)
%
%%

% Calculate bin information from given
r_pix = r / mm_per_pixel;                       % rail radius in pixels
[dims,rows,cols,p] = init_bin(frame{1},r_pix,angle);

% Initialize bin and mask
bin = zeros(dims(1),dims(2),'double');
mask = zeros(dims(1),dims(2),'like',bin);

% Create frame mask to track number of intersections at every pixel
frame_mask = ones(size(frame{1},1),size(frame{1},2),'like',bin);

for i = 1:length(frame)
    % Initialize bin and mask
    bin2add = zeros(dims(1),dims(2),'double');
    mask2add = bin2add;
    
    % Add frame to bin
    bin2add(rows(:),cols(:)) = bin2add(rows(:),cols(:)) + frame{i};
    
    % Plot B-scan before rotated
%     figure
%     subplot(1,2,1);
%     imagesc(bin2add); hold on
%     plot(p(2),p(1),'r+')
%     title('B-scan (Before Rotation)')
%     axis equal
    
    % Add mask frame to mask
    mask2add(rows(:),cols(:)) = mask2add(rows(:),cols(:)) + frame_mask;
    
    % Rotate frame around point in bin and mask
    bin2add = rotateAround(bin2add,p(1),p(2),-angle(i),method);
    mask2add = rotateAround(mask2add,p(1),p(2),-angle(i),method);
    
    % Plot rotated B-scan
%     subplot(1,2,2);
%     imagesc(bin2add); hold on
%     plot(p(2),p(1),'r+')
%     title('B-scan (Rotated)')
%     axis equal
    
    % Add to bin
    bin = bin + bin2add;
    mask = mask + mask2add;
end

% Average pixel values - pixel value divided by number of times B-scan
% intersected pixel
mask(mask == 0) = 1;
bin = abs(bin ./ mask);  % average pixel values


% Grayscale image
% subplot(1,2,1);
% imshow(uint8(bin)); axis on; hold on
% colorbar
% plot(p(2),p(1),'r+','MarkerSize',20);   % rotation point
% title('Reconstruction in Yaw (Grayscale)')
% xlabel('X (mm)')
% ylabel('Y (mm)')
% xticks(linspace(1,size(bin,2),10));
% xticklabels(linspace(1,size(bin,2) * mm_per_pixel,10));
% yticks(linspace(1,size(bin,1),10));
% yticklabels(linspace(1,size(bin,1) * mm_per_pixel,10));
% daspect([1 1 1])

% Colored image
% ax2 = subplot(1,2,2);
% ax2 = figure;
% imagesc(uint8(bin)); hold on
% colormap(ax2,parula)
% colorbar
% plot(p(2),p(1),'r+','MarkerSize',20);   % rotation point
% title('Reconstruction in Yaw (Colors)')
% xlabel('X (mm)')
% ylabel('Y (mm)')
% set(gca,'TickDir','out');
% xticks(linspace(0,size(bin,2),10));
% xticklabels(linspace(0,size(bin,2) * mm_per_pixel,10));
% xtickangle(45)
% yticks(linspace(0,size(bin,1),10));
% yticklabels(linspace(0,size(bin,1) * mm_per_pixel,10));
% daspect([1 1 1])
% axis on

end

%% Bin initialization
function [dims,rows,cols,p] = init_bin(frame,r,angle)
    % Dimensions of frame
    [fh,fw] = size(frame);
    
    % Find minimum bounding box
    center = floor([fw / 2, fh - r]);
    rows = [1 fh];
    cols = [1 fw];
    [bbox_x, bbox_y, center, x, y] = rot_bbox(cols, rows, center, angle);
    
    % Bin size (indexing to be 1:rows and 1:cols
    dims = [bbox_y(2), bbox_x(2)];   % [rows, cols]
    
    % Redefine points in the imshow coordinate axis - origin at top left
    % corner (positive x -> right, positive y -> down)
    
    % Center of rotation
    p = [dims(1) - center(2) + 1, center(1)];
    
    % Default frame in rows and cols
    rows = dims(1) - y(2) + 1:dims(1) - y(1) + 1;
    cols = x(1):x(2);
    
end
