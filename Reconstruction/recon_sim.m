clear; clc; close all

%% Reconstruction with simulated B-scans
% This script creates a set of binary images to simulate B-scans obtained
% using the design rail system. The images are accompanied with
% corresponding pose information. The reconstruction algorithm is used to
% reconstruct the frames into a volume and is compared to a ground truth
% volume.

%% PART 1: FAT SIDE ROTATION
% Create shape
shape = {'cube','sphere','cuboid'};
shapenum = 1;
[vol,center,corners] = createvol(shape{shapenum});

% Plot volume data
volIdx = find(vol);
[volX,volY,volZ] = ind2sub(size(vol),volIdx);
% scatter3(volX,volY,volZ,5,vol(volIdx),'o'); hold on

sz = size(vol);    % grid size
[X,Y,Z] = meshgrid(1:sz(1),1:sz(2),1:sz(3));
volX_re = reshape(volX,50,50,50);
volY_re = reshape(volY,50,50,50);
volZ_re = reshape(volZ,50,50,50);

% Slice location/position info
noSlices = 9;
angles = linspace(-45, 45, noSlices);
r = 20;
pose = {};

[ydef,zdef] = meshgrid(min(volY):max(volY), min(volZ):max(volZ));
planeSz = [max(volZ) - min(volZ) + 1, max(volY) - min(volY) + 1];
xdef = center(1) * ones(planeSz(1), planeSz(2));
method = {'linear','cubic','nearest'};

% Plot slices
subplot(1,2,1); hold on
frame = {};
    
for i = 1:noSlices
    [ydef,zdef] = meshgrid(min(volY):max(volY), min(volZ):max(volZ));
    
    % Rotation and translation offsets
    rot = zdef * sind(angles(i));
    dx = r * sind(abs(angles(i)));
    dz = r * (1 - cosd(abs(angles(i)))); 
    
    if angles(i) > 0
        xsurf = round(xdef + rot + dx);
    else
        xsurf = round(xdef + rot -dx);
    end
    
    zsurf = round(zdef - dz);
    
    % Remove invalid indices
    rm_idx = find(any(zsurf < 1, 2));
    xsurf(rm_idx,:) = [];
    ydef(rm_idx,:) = [];
    zsurf(rm_idx,:) = [];
    
    % Plot surfaces - differentiate between intersecting and
    % non-intersecting area
    surf_points = [xsurf(:),ydef(:),zsurf(:)];
    [intr,nonint] = findintersect(shape{shapenum},[volX,volY,volZ], ...
        surf_points);
    
    plane = zeros(size(ydef));
    plane(surf_points(intr,2)-min(surf_points(intr,2))+1,surf_points(intr,3)) = 255;
    plane = rot90(plane);
    
    scatter3(surf_points(intr,1),surf_points(intr,2),surf_points(intr,3),2,'b');
    scatter3(surf_points(nonint,1),surf_points(nonint,2),surf_points(nonint,3),2,'r');

    frame{end + 1} = plane;
end

% Plot cube/cuboid to show volume
if any(strcmpi(shape{shapenum},{'cube', 'cuboid'}))
    faces = [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8];
end

patch('Vertices',corners,'Faces',faces,'FaceColor','green','FaceAlpha',0.2);

title(['Slices of simulated ',shape{shapenum}])
legend('Intersecting','Non-Intersecting')
axis([0 sz(1) 0 sz(2) 0 sz(3)])
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
view(3)

% Plot reconstructed volume
bin = fillbin_thick(frame,angles,r,0);

figure; hold on
[x,y,z] = ind2sub(size(bin),find(bin));
scatter3(x,y,z,5,'filled')
title('Reconstruction of simulated frames')
axis([0 sz(1) 0 sz(2) 0 sz(3)])
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
view(3)

%% PART 2: THIN SIDE ROTATION
% Define ground truth cross section
box = zeros(51,51);
sq1 = ones(10,10);
sq1_center = [26,26];
noSlices = 9;
angles = linspace(-90, 90, noSlices);

frames = {};
for i = 1:noSlices
    frames{i} = 255 * ones(50);
end

% figure;
% bin_thin = fillbin_thin(frames,angle,100,'nearest');
% imshow(uint8(bin_thin))
% colormap jet
% colorbar

%% Sutherland-Hodgman Algorithm
% define the 6 planes of a box
% xrange = [-2 3];
% yrange = [-2 4];
% zrange = [-2 5];
% planes = [1  0  0  xrange(2);
%          -1  0  0 -xrange(1);
%           0  1  0  yrange(2);
%           0 -1  0 -yrange(1);
%           0  0  1  zrange(2);
%           0  0 -1 -zrange(1)];
% % intesect a bunch of random lines against the box
% for i=1:1000
%   pt1 = rand(1,3);
%   pt2 = 10*randn(1,3);
%   t = zeros(1,6);
%   for i=1:6
%     t(i) = Intersection(pt1,pt2,planes(i,:));
%   end
%   startpt = pt1;
%   endpt = pt1+min(t)*(pt2-pt1);
%   line([startpt(1) endpt(1)],[startpt(2) endpt(2)],[startpt(3) endpt(3)])
% end
% view(3)
% axis square
% 
% function t = Intersection(pt1,pt2,plane)
%   t = nan;
%   u = dot(pt1,plane(1:3)) - plane(4);    
%   v = dot(pt2,plane(1:3)) - plane(4);    
%   if (sign(u) ~= sign(v))
%       t = (0-u) / (v-u);
%   end
% end
