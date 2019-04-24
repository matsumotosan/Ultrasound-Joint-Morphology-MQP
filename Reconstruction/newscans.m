function [scans,shapebin] = newscans(frame_size,angles,r,shape,shape_size)
%NEWSCANS Create simulated US data
%
% NEWSCANS(FRAME_SIZE,ANGLES,RADIUS,SHAPE_SIZE,) - 

%%
% Frame dimensions
fh = frame_size(1); % height (no rows)
fw = frame_size(2); % width (no cols)

% Frame rotation dimensions
center = floor([fw / 2, fh - r]);   % center of rotation
rows = [1 fh];  % initial frame rows
cols = [1 fw];  % initial frame columns

% Find bounding box
[bbox_x, bbox_y, center, x, y] = rot_bbox(cols, rows, center, angles);

% Initialize bin
dims = [bbox_y(2), bbox_x(2)];
bin = zeros(dims(1),dims(2));    % initialize bin

% Default frame in rows and cols
rows = dims(1) - y(2) + 1:dims(1) - y(1) + 1;
cols = x(1):x(2);

% Center of rotation
p = [dims(1) - center(2) + 1, center(1)];

% Create second bin to contain shape
shapebin = zeros(dims(1),dims(2));

switch lower(shape)
    case 'square'
        shapebin(p(1) - shape_size / 2:p(1) + shape_size / 2, ...
                 p(2) - shape_size / 2:p(2) + shape_size / 2) = 1;
    case 'circle'
        [X,Y] = meshgrid(1:dims(2),1:dims(1));
        R = sqrt((X - p(2)) .^ 2 + (Y - p(1)) .^ 2);
        shapebin(R <= shape_size / 2) = 1;
    case 'composite'
        shapebin(p(1) - shape_size / 2:p(1) + shape_size / 2, ...
                 p(2) - shape_size / 2:p(2) + shape_size / 2) = 1;
        [X,Y] = meshgrid(1:dims(2),1:dims(1));
        dx = 10; dy = 10;
        R = sqrt((X - p(2) + dx) .^ 2 + (Y - p(1) + dy) .^ 2);
        shapebin(R <= shape_size / 2) = 1;
    otherwise
        error("SHAPE must be 'square' or 'circle'")
end

% Fill bin, rotate frame, overlap, rotate back, extract frame
scans = {};
% figure; hold on
for i = 1:length(angles)
    % Find intersection between framebin and shapebin
    framebin = bin;
    framebin(rows,cols) = framebin(rows,cols) + 2 * ones(fh,fw);
    framebin = rotateAround(framebin, p(1), p(2), angles(i));
    totalbin = framebin + shapebin;

%     subplot(1,length(angles),i);
%     imagesc(totalbin); hold on
%     plot(p(2),p(1),'r+')
    totalbin = rotateAround(totalbin, p(1), p(2), -angles(i));
    
    % Extract scan from totalbin
    [r,c] = find(totalbin > 1);
    s = totalbin(min(r):max(r),min(c):max(c));
    s(s < 2) = 2;
    s = s - 2;
    scans{i} = s(1:fh,1:fw);
end

% Show scans
% imagesc(scans{1}); hold on
% plot(p(2),p(1),'r+')


end

