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
[dims,rows,cols,p] = init_bin(frame{1},r_pix);  % bin info
bin = zeros(dims(1),dims(2),'double');          % initialize bin

for i = 1:length(frame)
    % Initialize bin
    bin2add = zeros(dims(1),dims(2),'double');
    
    % Add frame
    bin2add(rows(:),cols(:)) = bin2add(rows(:),cols(:)) + frame{i};
    
    % Rotate frame around point
    bin2add = rotateAround(bin2add,p(1),p(2),angle(i),method);
    
    % Add to bin
    bin = bin + bin2add;
end

% Average pixel values
% bin = floor(bin / length(frame));

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
ax2 = figure;
imagesc(uint8(bin)); hold on
colormap(ax2,parula)
colorbar
plot(p(2),p(1),'r+','MarkerSize',20);   % rotation point
title('Reconstruction in Yaw (Colors)')
xlabel('X (mm)')
ylabel('Y (mm)')
set(gca,'TickDir','out');
xticks(linspace(0,size(bin,2),10));
xticklabels(linspace(0,size(bin,2) * mm_per_pixel,10));
xtickangle(45)
yticks(linspace(0,size(bin,1),10));
yticklabels(linspace(0,size(bin,1) * mm_per_pixel,10));
daspect([1 1 1])
axis on

end

%%
function [dims,rows,cols,point] = init_bin(frame,r)
    % Dimension of frame
    [fh,fw] = size(frame);
    r_max = ceil(rssq([r, fw / 2]));
    
    % Reconstruction window size
    dims = ceil([r_max + fw / 2, 2 * r_max]);

    % Index of initial image insertion in placeholder bin
    rows = ceil(r_max - r + 1:r_max - r + fh);
    cols = ceil(r_max - fw / 2 + 1:r_max + fw / 2); 
    
    % Point to rotate around
    point = [2 * r_max - r,r_max];
end
