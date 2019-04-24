clear; clc; close all

%% Reconstruction with simulated B-scans
% This file contains scripts to demonstrate reconstruction in pitch and
% yaw. Reconstruction in pitch has not been tested, but reconstruction in
% yaw has been tested with some simple geometries. Run each section
% individually.

%% RECONSTRUCTION IN PITCH - DEMO
% This script creates a set of binary images to simulate B-scans obtained
% using the design rail system. The images are accompanied with
% corresponding pose information. The reconstruction algorithm is used to
% reconstruct the frames into a volume and is compared to a ground truth
% volume.

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
bin = fillbin_pitch(frame,angles,r,0);

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

%% RECONSTRUCTION IN YAW - BIN FILL DEMO
% % clear; clc; close all
% % 
% % % Rail and frame data
% % setting = load('us_setting_map.mat');
% % depth = 70;                             % US setting (mm)
% % win_data = setting.M(depth);            % load setting
% % mm_per_pixel = win_data(5);             % frame mm/pixel resolution
% % r = 70;                                 % rail radius (mm)
% % 
% % % Define slices
% % noSlices = 3;
% % angles = linspace(-45, 45, noSlices);
% % method = {'nearest','bilinear','bicubic'};
% % 
% % % Create fake frames
% % frames = {};
% % nrows = win_data(4) - win_data(3);
% % ncols = win_data(2) - win_data(1);
% % for i = 1:noSlices
% %     frames{i} = 100 * ones(nrows,ncols);
% % end
% % 
% % % Reconstruction in yaw
% % bin_thin = fillbin_yaw(frames,angles,r,mm_per_pixel,method{3});
%
% Compare interpolation methods
% for i = 1:length(method)
%     bin_thin = fillbin_thin(frames,angles,50,method{i});
%     ax = subplot(1,3,i); hold on
%     imshow(uint8(bin_thin)); axis on; hold on
%     colormap(ax,parula)
%     colorbar
%     caxis([0 max(max(bin_thin))])
%     title(method{i});
%     xlabel('Horizontal')
%     ylabel('Depth')
% end

%% RECONSTRUCTION IN YAW - SIMULATED B-SCAN
% Simulated B-scan reconstruction testing script. B-scans of simple
% geometries (square, circle, or a composite of the two) can be created and
% used for reconstruction in yaw. The expected image and reconstructed
% image are compared using two metrics: 2D cross correlation (normxcorr2)
% and the Structural Similarity Index (SSIM). The reconstruction algorithm
% (fillbin_yaw) requires the B-scans and their corresponding angles (among
% other arguments). Testing can be performed with ideal data. Random noise
% can be added artificially to the 'angles' input for the algorithm to
% simulate reconstruction with random errors that the IMU may experience
% during data acquisition.

clear; clc; close all

% Simulated B-scan parameters
frameSz = [60,20];
N = 20;
angles = linspace(-90,90,N);
radius = 30;
shape = {'square','circle','composite'};
shapeSz = 20;

% Introduce errors to angle
max_off = 0.96;
angles_err = angles + (-max_off + 2 * max_off * rand(1,N));

% Reconstruct
mm_per_pixel = 1;
method = {'Nearest','Bilinear','Bicubic'};

figure; hold on
for i = 1:length(method)
    % Simulated B-scans of shape
    [scans,shapebin] = newscans(frameSz,angles,radius,shape{2},shapeSz);
    
    % Reconstruct shape from simulated B-scans
    bin_yaw = fillbin_yaw(scans,angles,radius,mm_per_pixel,method{i});
    
    % Reconstruct shape from simulated B-scans (errors in angles)    
    bin_yaw_err = fillbin_yaw(scans,angles_err,radius,mm_per_pixel,method{i});
    
%     % Plot reconstructed shape and list SSIM and normxcorr2
%     ssimval = ssim(shapebin,bin_yaw);
%     xcorr_2d = max(max(normxcorr2(shapebin,bin_yaw)));
%
%     subplot(2,2,i+1); hold on
%     imagesc(bin_yaw_err)
%     title({[method{i}]; [' (SSIM:' num2str(round(ssimval,3)), ')']; ...
%           ['(normxcorr2:' ...
%           num2str(round(xcorr_2d,3)) ')']})
%     colorbar
%     axis equal tight
    
    % Compare reconstruction with random errors introduced
    [rows, cols] = size(shapebin);
    ssimval_err = ssim(shapebin,imresize(bin_yaw_err, [rows cols]));
    xcorr_2d_err = max(max(normxcorr2(shapebin,bin_yaw_err)));
    
    subplot(2,2,i+1)
    imagesc(bin_yaw_err)
    title({[method{i}]; [' (SSIM:' num2str(round(ssimval_err,3)), ')']; ...
          ['(normxcorr2:' ...
          num2str(round(xcorr_2d_err,3)) ')']})
    colorbar
    axis equal tight
end

% Plot original shape
subplot(2,2,1)
imagesc(shapebin)
title('Original')
colorbar
axis equal tight

% subplot(1,3,3)
% fz = imfuse(shapebin,bin_yaw,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
% imshow(fz)
% title('Fused')
% 
% % Image similarity
% figure
% C = normxcorr2(shapebin,bin_yaw);
% surf(C); shading flat
% title('Normalized Cross Correlation Coefficient')
% 
% err = immse(shapebin,bin_yaw);      % mean-squared error
% [peaksnr, snr] = psnr(shapebin, bin_yaw);  % peak signal-to-noise ratio


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
