function bin = fillbin_yaw(frame,angle,r,mm_per_pixel,method)
%FILLBIN_YAW - Distribute pixels to representative voxels
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

r_pix = r / mm_per_pixel;                       % rail radius in pixels
[dims,r,c,p] = init_bin(frame{1},r_pix,angle);  % bin info

bin = zeros(dims(1),dims(2),'double');          % initialize bin
mask = zeros(dims(1),dims(2),'like',bin);       % initialize mask
frame_mask = ones(size(frame{1},1),size(frame{1},2),'like',bin);

% Frame insertion position
rows = r(1):r(2);
cols = c(1):c(2);

for i = 1:length(frame)
    % Initialize bin and mask
    bin2add = zeros(dims(1),dims(2),'double');
    mask2add = bin2add;
    
    % Add frame to bin
    bin2add(rows(:),cols(:)) = bin2add(rows(:),cols(:)) + frame{i};
    
    figure
    imagesc(bin2add)
    
    % Add mask frame to mask
    mask2add(rows(:),cols(:)) = mask2add(rows(:),cols(:)) + frame_mask;
    
    % Rotate frame around point in bin and mask
    bin2add = rotateAround(bin2add,p(1),p(2),angle(i),method);
    mask2add = rotateAround(mask2add,p(1),p(2),angle(i),method);
    
    imagesc(bin2add)
    
%     figure
%     imagesc(bin2add)
    
    % Add to bin
    bin = bin + bin2add;
    mask = mask + mask2add;
    
    imagesc(bin)
end

% Average pixel values
mask(mask < 1) = 1;
bin = bin ./ mask;

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

%%
function [dims,rows,cols,center] = init_bin(frame,r,angle)
    % Dimensions of frame
    [fh,fw] = size(frame);
    
    % Find minimum bounding box
    center = floor([fw / 2, r - fh]);
    rows = [1 fh];
    cols = [1 fw];
    [bbox_x, bbox_y] = rot_bbox(cols, rows, center, angle, 'y');
    
    % Calculate offset to convert all coordinates relative to global origin
    shift_x = -min(bbox_x) + 1;
    if min(bbox_y) < center(2)
        shift_y = -min(bbox_y) + 1;
    else
        shift_y = -center(2) + 1;
    end
    
    % Shift center
    center = center + [shift_x, shift_y];
    
    % Shift default rows and columns of frame
    rows = rows + shift_y;
    cols = cols + shift_x;
    
    % Align bottom left corner of bbox with (0,0)
    bbox_x = bbox_x + shift_x;
    bbox_y = bbox_y + shift_y;
    dims = [max(bbox_y) max(bbox_x)];
    
%     % Plot newly calculated center, default frame location, and bbox
%     % Shifted center of rotation
%     figure; hold on 
%     plot(center(1), center(2), 'r+');
%     
%     % Default frame location
%     k = convhull(x,y);
%     plot(x(k),y(k),'-ro');
%     
%     % Shifted minimum bounding box
%     
%     k = convhull(bbox_x, bbox_y);
%     plot(bbox_x(k), bbox_y(k), '-k*');
%     
%     % All rotated frames
%     for i = 1:length(angle)
%         k = convhull(x2(i,:), y2(i,:));
%         plot(x2(i,k),y2(i,k),'--bo');
%     end
% 
%     title('Minimum bounding box of frames (shifted)');
%     xlabel('X'); ylabel('Y');
%     legend('Center of Rotation', 'Original', 'Bounding Box', 'Rotated Frames');
%     axis equal
    
%     if (r > fh)
%         r_max = ceil(rssq([r, fw / 2]));
%     else
%         r_max = ceil(rssq([fh, fw / 2]));
%     end
%     
%     % Reconstruction window size
%     dims = ceil([r_max + fw / 2, 2 * r_max]);
% 
%     % Index of initial image insertion in placeholder bin
%     rows = ceil(r_max - r + 1:r_max - r + fh);
%     cols = ceil(r_max - fw / 2 + 1:r_max + fw / 2); 
%     
%     % Point to rotate around
%     point = [2 * r_max - r,r_max];
end
