function [rectx, recty] = rot_bbox(x,y,center,theta,varargin)
%ROT_BBOX - Calculate corners of minimum bounding box of rectangle
%rotated by theta
%
%   ROT_BBOX(X,Y,CENTER,THETA) - returns coordinates of bounding box of
%   points four sets of points (X,Y) rotated by THETA (degrees) around 
%   CENTER (X,Y). Default is to not plot result. RECTX and RECTY are both 
%   1-by-2 vectors.
%
%   ROT_BBOX(X,Y,CENTER,THETA,OPT) - if OPT is 'y', plot original rectangle
%   minimum bounding box, and all rotated rectangles
%
%%

% Process input args
plotting = 0;
if nargin > 4
    if strcmpi(varargin{1}, 'y')
        plotting = 1;
    end
end

% Preprocess inputs
x = [x(1) x(2) x(1) x(2)];
y = [y(1) y(1) y(2) y(2)];
theta = theta(:);

% Calculate x and y for all rotation by theta around center
x2 = center(1) + (x - center(1)) .* cosd(theta) + (y - center(2)) .* sind(theta);
y2 = center(2) - (x - center(1)) .* sind(theta) + (y - center(2)) .* cosd(theta);

% Bounding box x and y - if not integer take ceiling/floor to contain all
% points
rectx = [floor(min(min(x2))) ceil(max(max(x2)))];   % [xmin xmax]
recty = [floor(min(min(y2))) ceil(max(max(y2)))];   % [ymin ymax]

if (plotting)
    % Plot center of rotation
    figure; hold on 
    plot(center(1), center(2), 'r+');
    
    % Plot original rectangle
    k = convhull(x,y);
    plot(x(k),y(k),'-ro');
    
    % Plot minimum bounding box
    bbox_x = [rectx(1) rectx(1) rectx(2) rectx(2)];
    bbox_y = [recty(1) recty(2) recty(1) recty(2)];
    k = convhull(bbox_x, bbox_y);
    plot(bbox_x(k), bbox_y(k), '-k*');
    
    % Plot rotated rectangles
    for i = 1:length(theta)
        k = convhull(x2(i,:), y2(i,:));
        plot(x2(i,k),y2(i,k),'--bo');
    end

    title('Minimum bounding box of frames');
    xlabel('X'); ylabel('Y');
    legend('Center of Rotation', 'Original', 'Bounding Box', 'Rotated Frames');
    axis equal
end

end