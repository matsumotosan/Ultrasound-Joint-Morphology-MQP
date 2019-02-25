%% Create colored video with input file names option

close all; clear all; clc;
name_of_mat = input('Name of .mat file?: ','s');
name_of_mat = sprintf('%s.mat',name_of_mat);
load(name_of_mat)

name_of_vid = input('Name of .avi file?: ','s');
name_of_vid = sprintf('%s.avi',name_of_vid);

FrameRate = input('Frame Rate?: ');

writerObj = VideoWriter(name_of_vid);
writerObj.FrameRate = FrameRate;

% open the video writer
open(writerObj);

% write the frames to the video
for u=1:length(frames)
    % convert the image to a frame
    frame = im2frame(frames{u});
    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);


%% convert to grayscale attempt
close all; clear all; clc;
name_of_mat = input('Name of .mat file?: ','s');
name_of_mat = sprintf('%s.mat',name_of_mat);
load(name_of_mat)

grayFrames = cell(1,length(frames));
for i = 1:length(frames)
    gray = frames{i};
    grayFrames{i} = rgb2gray(gray);    
end

name_of_vid = input('Name of .avi file?: ','s');
name_of_vid = sprintf('%s.avi',name_of_vid);

FrameRate = input('Frame Rate?: ');

writerObj = VideoWriter(name_of_vid);
writerObj.FrameRate = FrameRate;

% open the video writer
open(writerObj);

% write the frames to the video
for u=1:length(frames)
    % convert the image to a frame
    frame = im2frame(frames{u});
    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);

 vid = VideoReader(name_of_vid);
 numImgs = get(vid, 'NumberOfFrames');
 vid_frames = read(vid);
 
 grayvidname = sprintf('gray_%s',name_of_vid);
 
 obj=VideoWriter(grayvidname);
 open(obj);

 for i=1:numImgs
     movie(i).cdata=rgb2gray(vid_frames(:,:,:,i));
     movie(i).colormap=gray;
 end

 writeVideo(obj,movie);
 close(obj);

%% Gif Attempt


h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for n = 1:0.5:5
    % Draw plot for y = x.^n
    x = 0:0.01:1;
    y = x.^n;
    plot(x,y) 
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
  end

