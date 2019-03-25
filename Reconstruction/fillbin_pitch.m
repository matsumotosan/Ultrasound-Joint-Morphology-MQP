function bin = fillbin_pitch(frame,pose,r,pad)
%FILLBIN_THICK - Distribute pixels to representative voxels
%
%   BIN = FILLBIN(FRAME,POSE,R,PAD) rotates and translates every from
%   from FRAME according to the corresponding POSE and R and fills BIN. The
%   OPT specifies the orientation the frame is rotated. PAd can be specified 
%   to add thickness to each of the frames to fill more of BIN.    
%
%   BIN = FILLBIN(FRAME,POSE,R,T)
%
%   BIN = FILLBIN(___,METHOD)
%

if mod(size(frame{1},2),2) == 0
    for i = 1:length(frame)
        frame{i} = padarray(frame{i},1,0,'pre');
    end
end

axis = [1 0 0];
[dims,def] = init_bin(frame{1},pose,r,pad);     % calculate bin size
bin = zeros(dims(1),dims(2),dims(3),'double');  % initialize bin

% Add each frame to bin
for i = 1:length(frame)
    % Initialize bin
    vox2add = zeros(dims(1),dims(2),dims(3),'double');
    
    % Add frame
    vox2add(def{1},def{2},def{3}) = squeeze(vox2add(def{1},def{2},def{3})) + flip(frame{i},1)';
    
    figure; hold on
    [x,y,z] = ind2sub(size(vox2add),find(vox2add));
    scatter3(x,y,z,5,'filled')

    % Rotate frame at center
    vox2add = imrotate3(vox2add,pose(i),axis,'nearest','crop');

    [x,y,z] = ind2sub(size(vox2add),find(vox2add));
    scatter3(x,y,z,5,'filled')
    
    % Translate frame
    dx = r * sind(pose(i));
    dz = r * (1 - cosd(pose(i)));
    
    vox2add = imtranslate(vox2add, [0 dx dz],'nearest');
     
    [x,y,z] = ind2sub(size(vox2add),find(vox2add));
    scatter3(x,y,z,5,'filled')
    title(['Frame ', num2str(i), ' Transformation Process'])
    xlabel('x'); ylabel('y'); zlabel('z');
    legend('Default','Rotated','Translated')
    
    % Add to bin
    bin = bin + vox2add;
end

%%
function [dims,def] = init_bin(frame,pose,r,pad)

    % Pixel dimensions of frames
    [fh,fw] = size(frame);

    % Number of rows (x-dim)
    rows = ceil((r + fh) * (sind(abs(max(pose))) + sind(abs(min(pose)))));

    % Height (z-dim)
%     zmax = r * sind(90 - min(abs(pose)));
%     zmin = r * sind(90 - max(abs(pose)));
    height = ceil(rows + fw / 2);
    
    % Dimension of bounding box
    dims = [rows,fw,height];

    % Indices for slice insertion in placeholder bin
    def = {round(rows / 2) - pad:round(rows / 2) + pad, 1:fw, ...
        round(height / 2) - (fh - 1) / 2:round(height / 2) + (fh - 1) / 2};
end

end

