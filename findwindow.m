%% Find cropping window to remove non-image related pixels from scans
% Load frames
addpath(genpath('/Users/Shion/Dropbox/MQP-US Media'))
file = '5_cm.mat';    % Shion
im = load(file);
im = rgb2gray(im.frames{1});
imbw = imbinarize(im,'adaptive');

% Show original and binarized frames
figure; hold on
subplot(1,2,1)
imshow(im)
title('Original')
subplot(1,2,2)
imshow(imbw)
title('Binarized');

% Find coordinates of cropping window
B = bwboundaries(imbw);
imshow(imbw);hold on
visboundaries(B)
axis on

% Segment using active contours (snake) algorithm
% Specify initial contour location close to object to be segmented
% mask = false(size(imbw));
% mask(60:420,270:350) = true;
% figure
% imshow(imbw); hold on; axis on;
% visboundaries(mask,'Color','b')
% bd = activecontour(imbw,mask,50,'edge');
% visboundaries(bd,'Color','r');
% title('Blue - Initial Contour, Red - Final Contour')

