function bin = fillbin(frame,pose,r,scale,opt,pad)
%FILLBIN - distribute pixels to representative voxels
%
%   BIN = FILLBIN(FRAME,POSE,R,T) rotates and translates every from
%   from FRAME according to the corresponding POSE and R and fills BIN. The
%   OPT specifies the orientation the frame is rotated. T can be specified 
%   to add thickness to each of the frames to fill more of BIN.    
%
%   BIN = FILLBIN(FRAME,POSE,R,T)
%
%   BIN = FILLBIN(___,METHOD)
%
%
%%
% [opt,t] = parseinputs(varargin{:});

if strcmp(opt,'yaw')
    axis = [1 0 0];
elseif strcmp(opt,'pitch')
    axis = [0 1 0];
else
    error('opt must be either yaw or pitch')
end

[dims,def] = init_bin(frame{1},pose,r,scale,pad);   % calculate bin size
bin = zeros(dims(1),dims(2),dims(3),'uint8');       % initialize bin

% Add each frame to bin
for i = 1:length(frame)
    % Initialize bin
    vox2add = zeros(voxSz(1),voxSz(2),voxSz(3),'uint8');
    
    % Add frame
    %vox2add(defX-1:defX+1,:,defZ) = squeeze(vox2add(defX-1:defX+1,:,defZ)) + repmat(frame{i}',[1 1 3]);
%     vox2add(defX,:,defZ) = vox2add(defX,:,defZ) + shiftdim(repmat(frame{i}',[1 1 1]),2);
vox2add(def(1),:,def(2)) = squeeze(vox2add(def(1),:,def(2))) + frame{i}';
    
    % Rotate frame
    vox2add = imrotate3(vox2add,pose(i),axis,'nearest','crop');
%     vox2add = imrotate3_fast(vox2add,[pose(i) 0 0],'nearest');

    % Translate frame
    dx = r * sind(abs(pose(i)));
    dz = r * (1 - cosd(abs(pose(i))));
    
    if pose(i) < 0
        vox2add = imtranslate(vox2add, [0 -dx -dz]);
    elseif pose(i) > 0
        vox2add = imtranslate(vox2add, [0 dx -dz]);
    end
    
    % Add to bin
    bin = bin + vox2add;
end

% Find bounding box
% [x,y,z] = ind2sub(size(bin),find(bin));
% [~,corners,~,~,~] = minboundbox(x,y,z);
% 
% corners = round(corners);
% 
% xrange = min(corners(:,1)):max(corners(:,1));
% yrange = min(corners(:,2)):max(corners(:,2));
% zrange = min(corners(:,3)):max(corners(:,3));
% bin = bin(xrange,yrange,zrange);

% Average all voxels
% bin = bin / length(frame);

% Plot result
% idx = find(bin);
% [a,b,c] = ind2sub(size(bin),idx);
% figure
% scatter3(a,b,c); hold on
% 
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% 
% xlim([0 voxSz(1)])
% ylim([0 voxSz(2)])
% zlim([0 voxSz(3)])

% function [t,opt] = parseinputs(varargin)
%     
%     % Check inputs
%     narginchk(4,6);
%     
% end
    
    function [dims,def] = init_bin(frame,pose,r,scale,pad)
        
        % Pixel dimensions of frames
        [fh,fw] = size(frame);
        
        % Range of x in real dimensions
        xrange = (r + fh / 2) * (sin(abs(max(pose))) + sin(abs(min(pose))));
        
        if min(abs(pose))
            zmax = (r + fh / 2) * cos(min(abs(pose(pose ~= 0))));
        else
            zmax = r + fh / 2;
        end
        
        zmin = (r - fh / 2) * cos(max(abs(pose)));
        
        % Range of z in real dimensions
        zrange = zmax - zmin;
        
        % No yrange because y-dim is width of frame
        dims = ceil([xrange / scale, fw, zrange / scale]);
        
        % Indices for slice insertion in placeholder
        def = {round(dims(1) / 2) - pad:round(dims(1) / 2) + pad, ...
               dims(3) - fh / scale:dims(3)};
    end

end

