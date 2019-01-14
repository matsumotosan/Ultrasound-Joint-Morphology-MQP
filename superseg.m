function video = superseg(data, n, comp)
%SUPERSEG Perform superpixel oversegmentation on every frame in video
%   Reads every frame in video and performs superpixel oversegmentation.
%   Once complete, assign average binary value of of superpixel to
%   corresponding superpixel.
%   
%   Input:  file - file name of video
%   
%   Output: video - cell array of segmented frames
%
%   Video/Image processing toolbox needed

v = VideoReader(data);
cropreg = [30 87 515 250];  % cropping region for entire video
comp = 5;                   % compactness of superpixellation
video = {};                 % initialize video

while hasFrame(v)
    f = readFrame(v);                   % read next frame of video
    f = rgb2gray(imcrop(f, cropreg));   % crop and convert to grayscale
    f = imadjust(f,[0.01,0.1],[]);      % increase contrast of frame
    
    count = 1;
    
    % Superpixel oversegmentation
    for i = n
        [L,N] = superpixels(f,i,'Compactness',comp(count));

        f_out = zeros(size(f),'like',f);    % initialize zeros for output image
        idx = label2idx(L);                 % linear indices of superpixel regions

        % Assign average of superpixel to corresponding superpixel
        for labelVal = 1:N
            grayIdx = idx{labelVal};
            f_out(grayIdx) = mean(f(grayIdx));
        end
        
        f = f_out;
%         imshow(f)
%         pause(0.1);
    end
    
    % Compare image before/after segmentation
%     subplot(1,2,1)
%     imshow(f)
%     subplot(1,2,2)
%     imshow(f_out);
%     pause(0.1)
    
%     imshow(f)
%     pause(0.1)

    % Save frame
    video{end+1} = f;
end

end

