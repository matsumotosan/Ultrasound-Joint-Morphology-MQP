clear; clc; close all

%% Reconstruction with simulated B-scans
% This script creates a set of binary images to simulate B-scans obtained
% using the design rail system. The images are accompanied with
% corresponding pose information. The reconstruction algorithm is used to
% reconstruct the frames into a volume and is compared to a ground truth
% volume.

%% PART 1: NO ERRORS IN IMAGE OR POSE
% Create shape
shape = {'cube','sphere','cuboid'};
shapenum = 1;
[vol,center] = createvol(shape{shapenum});

% Plot volume data
volIdx = find(vol);
[volX,volY,volZ] = ind2sub(size(vol),volIdx);
scatter3(volX,volY,volZ,5,vol(volIdx),'o'); hold on

sz = size(vol);    % grid size
[X,Y,Z] = meshgrid(1:sz(1),1:sz(2),1:sz(3));

% Slice location/position info
% fh = 30;
% fw = 50;
noSlices = 1;
% angles = linspace(-90, 15, noSlices);
angles = 45;
r = 20;
pose = {};

[ysurf,zsurf] = meshgrid(min(volY):max(volY), min(volZ):max(volZ));
planeSz = [max(volZ) - min(volZ) + 1, max(volY) - min(volY) + 1];
xdef = center(1) * ones(planeSz(1), planeSz(2));
method = {'linear','cubic','nearest'};

% Plot slices
% figure; hold on
lgd = {'Volume'};
% shp = alphaShape(volX,volY,volZ);
% plot(shp,'FaceAlpha',0.4)

for i = 1:noSlices
    % Rotate and translate
    offset = zsurf * sin(angles(i));
    xsurf = round(xdef + offset);
    
    % Remove invalid indices
    rm_idx = find(isequal();
    
    % Plot surfaces - differentiate between intersecting and
    % non-intersecting area
    surf(xsurf,ysurf,zsurf)
%     slice(X,Y,Z,double(vol),xsurf,ysurf,zsurf);
    lgd{end + 1} = num2str(angles(i));
end

title(['Slices of simulated ',shape{shapenum}])
legend(lgd)
axis([0 sz(1) 0 sz(2) 0 sz(3)])
xlabel('X')
ylabel('Y')
zlabel('Z')

% Reconstruct
% bin = fillbin(s,pose,angles,10,1,'yaw',1);


%% PART 2: 
% [X,Y,Z] = meshgrid(-5:0.2:5);
% V = X.*exp(-X.^2-Y.^2-Z.^2);
% 
% [xsurf,ysurf] = meshgrid(-2:0.2:2);
% zsurf = xsurf.^2-ysurf.^2;
% slice(X,Y,Z,V,xsurf,ysurf,zsurf)


