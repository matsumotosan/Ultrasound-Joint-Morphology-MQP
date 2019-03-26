function [bbox_x, bbox_y, center, xnew, ynew] = rot_bbox(x,y,center,theta,varargin)
%ROT_BBOX - Calculate corners of minimum bounding box of rectangle
%rotated by theta around a point
%
%   ROT_BBOX(X,Y,CENTER,THETA) - returns coordinates of bounding box of
%   points four sets of points (X,Y) rotated by THETA (degrees) around
%   CENTER. Also returns coordinates of center of rotation, CENTER, and
%   default frame location for input frames XNEW and YNEW. All coordinates
%   are shifted such that the left corner of the bounding box is at (1,1)
%   on the x-y axis. Counterclockwise rotation defined as positive.
%
%   ROT_BBOX(X,Y,CENTER,THETA,OPT) - if OPT is 'y', plot original rectangle
%   minimum bounding box, center of rotation, and all rotated rectangles.
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
X = [x(1) x(2) x(1) x(2)];
Y = [y(1) y(1) y(2) y(2)];
theta = theta(:);

% Calculate x and y after rotation by theta around center
x2 = center(1) + (X - center(1)) .* cosd(-theta) + (Y - center(2)) .* sind(-theta);
y2 = center(2) - (X - center(1)) .* sind(-theta) + (Y - center(2)) .* cosd(-theta);

% Bounding box x and y - if not integer take ceiling/floor to contain all
% points
bbox_x = [floor(min(x2(:))), ceil(max(x2(:)))];
bbox_y = [floor(min([y(:); y2(:)])), ceil(max([y(:); y2(:)]))];

% Plot unshifted bounding box, center, default frame location
if (plotting)
    subplot(1,2,1); hold on
    plot_bbox
    title('Minimum bounding box of frames (Unshifted)');
    hold off
end

% Calculate offset to convert all coordinates relative to global origin
shift_x = -min(bbox_x) + 1;
if min(bbox_y) < center(2)
    shift_y = -min(bbox_y) + 1;
else
    shift_y = -center(2) + 1;
end

% Shift center
center = center + [shift_x, shift_y];

% Shift x and y of frame
X = X + shift_x;
Y = Y + shift_y;
xnew = unique(X);
ynew = unique(Y);

% Shift rotated frames
x2 = x2 + shift_x;
y2 = y2 + shift_y;

% Align bottom left corner of bbox with (0,0)
bbox_x = bbox_x + shift_x;
bbox_y = bbox_y + shift_y;

% Plot shifted bounding box
if (plotting)
    subplot(1,2,2); hold on
    plot_bbox
    title('Minimum bounding box of frames (Shifted)');
    hold off
end


%% Bounding box plotting
function plot_bbox
    plot(center(1), center(2), 'r+');

    % Plot original rectangle
    k = convhull(X,Y);
    plot(X(k),Y(k),'-ro');

    % Plot minimum bounding box
    rect_x = [bbox_x(1) bbox_x(1) bbox_x(2) bbox_x(2)];
    rect_y = [bbox_y(1) bbox_y(2) bbox_y(1) bbox_y(2)];
    k = convhull(rect_x, rect_y);
    plot(rect_x(k), rect_y(k), '-k*');

    % Plot all rotated rectangles
    for i = 1:length(theta)
        k = convhull(x2(i,:), y2(i,:));
        plot(x2(i,k),y2(i,k),'--bo');
    end

    xlabel('X'); ylabel('Y');
    legend('Center of Rotation', 'Default Frame', 'Bounding Box', 'Rotated Frames');
    axis tight
end

end