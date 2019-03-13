function [scans,shapebin,mv] = newscans(frame_size,angles,r,shape,shape_size)
%NEWSCANS Create simulated US data
%
% NEWSCANS(FRAME_SIZE,ANGLES,RADIUS,SHAPE_SIZE,) - 

%%
% Initialize bin
fh = frame_size(1);
fw = frame_size(2);

r_max = ceil(norm([r, fw / 2]));            % traversing radius
dims = ceil([r_max + fw, 2 * r_max]);       % bin size
bin = zeros(dims(1),dims(2));               % initialize bin

rows = ceil(r_max - r + 1:r_max - r + fh);          % default frame rows
cols = ceil(r_max - fw / 2 + 1:r_max + fw / 2);     % default frame columns

p = [2 * r_max - r,r_max];  % point to rotate around, also center of shape

% Create second bin with shape
shapebin = zeros(dims(1),dims(2));

switch lower(shape)
    case 'square'
        shapebin(p(1) - shape_size / 2:p(1) + shape_size / 2, ...
                 p(2) - shape_size / 2:p(2) + shape_size / 2) = 1;
    case 'circle'
        [X,Y] = meshgrid(1:dims(2),1:dims(1));
        R = sqrt((X - p(2)) .^ 2 + (Y - p(1)) .^ 2);
        shapebin(R <= shape_size) = 1;
    otherwise
        error("SHAPE must be 'square' or 'circle'")
end

% Fill bin, rotate frame, overlap, rotate back, extract frame
scans = {};
for i = 1:length(angles)
    framebin = bin;
    framebin(rows,cols) = framebin(rows,cols) + 2 * ones(fh,fw);
    framebin = rotateAround(framebin, p(1), p(2), angles(i));
    totalbin = framebin + shapebin;
%     figure
    imagesc(totalbin)
    totalbin = rotateAround(totalbin, p(1), p(2), -angles(i));
%     figure
%     imagesc(totalbin)
    mv(i) = getframe;
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

