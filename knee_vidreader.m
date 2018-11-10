clear; clc; close all;

% Trial script to perform offline image-processing
% Ensure video is of compatible file format (.mp4, .avi)
% Edge detection using MATLAB defined algorithms
% 
% Requires following toolboxes:
% -Image processing
% -Video processing
% -Computer vision toolbox
%% Read video
file = 'Analog_20181109_1531.mp4';

% Play video
videoFReader = vision.VideoFileReader(file);
videoPlayer = vision.VideoPlayer;
while ~isDone(videoFReader)
  videoFrame = videoFReader();
  videoPlayer(videoFrame);
  pause(0.1)
end

%% Edge detection
% Read specific frames
v = VideoReader(file);
frames = 1:10;
im = rgb2gray(read(v,4));           % convert image to grayscale
crop_rect = [29.5 86.5 513 370];    % cropping parameters
im = imcrop(im, crop_rect);
imshow(im);

% Perform edge detection with various algorithms
edges_canny = edge(im,'Canny');
edges_prewitt = edge(im,'Prewitt');
edges_roberts = edge(im,'Roberts');
edges_sobel = edge(im,'Sobel');

% Plot edge detection with basic methods
figure;
subplot(2,2,1);
imshow(edges_canny);
title('Canny');
subplot(2,2,2);
imshow(edges_prewitt);
title('Prewitt');
subplot(2,2,3);
imshow(edges_roberts);
title('Roberts');
subplot(2,2,4);
imshow(edges_sobel);
title('Sobel');

% Prewitt algorithm experiment
edges_prewitt_h = edge(im,'Prewitt','horizontal');

figure;
subplot(1,2,1);
imshow(edges_prewitt);
title('Prewitt');
subplot(1,2,2);
imshow(edges_prewitt_h);
title('Prewitt (horizontal)');

% Roberts algorithm experiment
edges_roberts_h = edge(im,'Roberts','horizontal');

figure;
subplot(1,2,1);
imshow(edges_roberts);
title('Roberts');
subplot(1,2,2);
imshow(edges_roberts_h);
title('Roberts (horizontal)');

% Sobel algorithm experiment
edges_sobel_h = edge(im,'Sobel','horizontal');

figure;
subplot(1,2,1);
imshow(edges_sobel);
title('Sobel');
subplot(1,2,2);
imshow(edges_sobel_h);
title('Sobel (horizontal)');
