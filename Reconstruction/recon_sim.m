clear; clc; close all

%% Reconstruction with simulated B-scans
% This script creates a set of binary images to simulate B-scans obtained
% using the design rail system. The images are accompanied with
% corresponding pose information. The reconstruction algorithm is used to
% reconstruct the frames into a volume and is compared to a ground truth
% volume.

%% PART 1: NO ERRORS IN IMAGE OR POSE
shape = {'cube','sphere','cuboid'};
shapenum = 1;
vol_true = createvol(shape{shapenum});

sz = size(vol_true);
[X,Y,Z] = meshgrid(1:sz(1),1:sz(2),1:sz(3));

pose = {};
figure; hold on
for i = 1:length(noSlices)
    slice(X,Y,Z,);
    poses{end + 1} = num2str(angles(i));
end

title(['Slices of simulated ',shape{shapenum}])
legend(poses)


%% PART 2: 