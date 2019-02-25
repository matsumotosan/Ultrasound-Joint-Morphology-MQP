function dims = rotminbox(corners,angles)
%ROTMINBOX Find minimum bounding box for rectangle with CORNERS rotated
%around origin
%
%   ROTMINBOX(CORNERS,ANGLES,POINT) finds the dimensions DIMS of the 
%   minimum bounding box containing a rectangle with CORNERS rotated at 
%   ANGLES around POINT. 

% Parameter check
if length(corners) ~= 4 || corners(1) > corners(2) || corners(3) > corners(4)
    error('CORNERS must be in the format [x1 x2 y1 y2])')
end

% Gather initial coordinates and rotation data
topleft = [corners(3) corners(1)];
topright = [corners(3) corners(2)];
ccw = 0:max(angles);
cw = 0:min(angles);

% Find corners of bounding box
xmin = floor(min());
xmax = ceil(max());

end

