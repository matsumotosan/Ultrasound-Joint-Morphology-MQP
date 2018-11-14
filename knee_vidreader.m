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
% videoFReader = vision.VideoFileReader(file);
% videoPlayer = vision.VideoPlayer;
% while ~isDone(videoFReader)d
%   videoFrame = videoFReader();
%   videoPlayer(videoFrame);
%   pause(0.1)
% end

%% Preprocessing
v = VideoReader(file);
frame = 10;                         % frames of interest
im = rgb2gray(read(v,frame));       % convert image to grayscale
crop_rect = [29.5 86.5 513 250];    % cropping parameters
im = imcrop(im, crop_rect);         % crop image
im_adjust = imadjust(im,[0.01 0.1],[]); % magnify contrast
imshow(im_adjust)

%% 2-D superpixel oversegmentation of images
n = 20;         % number of superpixels to create
comp = 15;      % compactness: higher value makes superpixels more regularly shaped
[L,N] = superpixels(im_adjust,n,'Compactness',comp);
mask = boundarymask(L);

% Compare to original image
subplot(1,2,1)
imshow(im_adjust)
subplot(1,2,2)
imshow(imoverlay(im_adjust,mask,'cyan'),'InitialMagnification',67);
title({strcat('Number of superpixels: ', num2str(n)), ...
    strcat('Compactness: ',num2str(comp))});

% Superpixel posterization
outputImage = zeros(size(im_adjust),'like',im_adjust);  % initialize zeros of similar type
idx = label2idx(L); % create cell array of linear indices for superpixel positions
numRows = size(im_adjust,1);
numCols = size(im_adjust,2);
for labelVal = 1:N
    grayIdx = idx{labelVal};
    outputImage(grayIdx) = mean(im_adjust(grayIdx));
end    

figure
imshow(outputImage,'InitialMagnification',67)
title('Superpixel posterization');

%% Function: superseg
% Perform superpixel oversegmentation on all frames of video
v = superseg(file);

% Play video
figure
for i = 1:length(v)
    imshow(v{i});
    pause(0.2);
end

%% Edge detection experiment
% Compare original cropped and adjusted crop
figure;
subplot(1,2,1)
imshow(im);
title('Raw Image');
subplot(1,2,2)
imshow(im_adjust)
title('Contrast Adjusted Image');

% Perform edge detection with various algorithms
edges_canny = edge(im_adjust,'Canny');
edges_prewitt = edge(im_adjust,'Prewitt');
edges_roberts = edge(im_adjust,'Roberts');
edges_sobel = edge(im_adjust,'Sobel');

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
edges_prewitt_h = edge(im_adjust,'Prewitt','horizontal');

figure;
subplot(1,2,1);
imshow(edges_prewitt);
title('Prewitt');
subplot(1,2,2);
imshow(edges_prewitt_h);
title('Prewitt (horizontal)');

% Roberts algorithm experiment
edges_roberts_h = edge(im_adjust,'Roberts','horizontal');

figure;
subplot(1,2,1);
imshow(edges_roberts);
title('Roberts');
subplot(1,2,2);
imshow(edges_roberts_h);
title('Roberts (horizontal)');

% Sobel algorithm experiment
edges_sobel_h = edge(im_adjust,'Sobel','horizontal');

figure;
subplot(1,2,1);
imshow(edges_sobel);
title('Sobel');
subplot(1,2,2);
imshow(edges_sobel_h);
title('Sobel (horizontal)');

%% Cavity fill experiment
im_fillholes = imfill(im_adjust, 'holes');
im_fill = imfill(im_adjust);
subplot(1,2,1)
imshow(im_fillholes);
title('imfill - holes');
subplot(1,2,2)
imshow(im_fill)
title('imfill - no input');

%% Image filtering - 2D Gaussian filter
% Gaussian image filtering - commonly used to reduce noise
sigma = 2;    % default: 0.5
im_gaussfilter = imgaussfilt(im_adjust,sigma);

% Compare to original contrast-enhanced image
subplot(1,2,1)
imshow(im_adjust)
title('Original - contrast enhanced');
subplot(1,2,2)
imshow(im_gaussfilter);
title(strcat('Gaussian filter (\sigma = ',num2str(sigma),')'));

% Try edge detection

%% Image filtering - bilateral filtering with Gaussian kernels


%% Hough transform experiment
[H,T,R] = hough(im_adjust);
imshow(imadjust(rescale(H)),'XData',T,'YData',R,...
      'InitialMagnification','fit');
title('Hough transform');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(gca,hot);

% Find peaks in Hough transform of image
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');

% Find lines and plot them
lines = houghlines(im_adjust,T,R,P,'FillGap',5,'MinLength',7);
figure, imshow(im_adjust), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if (len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
