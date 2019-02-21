function vol = createvol(varargin)
%CREATEVOLDATA Create 3D matrix to simulate volumetric data
%
%   VOL = CREATEVOL(SHAPE) creates a 3-D uint8 matrix in the SHAPE
%   specified in the center of a 100-by-100-by-100 matrix. If cube, side
%   length is 50. If cuboid side lengths are 20, 40, 60 in x-, y-, and
%   z-directions, respectively. If sphere, radius is set to 25. If
%   ellipsoid, radii in x-, y-, and z-directions are set to 10, 20, and 30,
%   respectively.
%   
%   VOL = CREATEVOL(SHAPE,DIM,LIM) creates a 3-D uint8 matrix in the SHAPE
%   specified with dimensions DIM in a box of size LIM. If cube, specify 
%   side length. If cuboid, specify side lengths as [sx sy sz]. If sphere,
%   specify radius. If ellipsoid, specify radii as [rx ry rz]. 
%
%   VOL = CREATEVOL(___,GRADIENT) creates volume with color gradient in
%   direction specified by GRADIENT
%   
%
%%
[D,lim,grad] = parseinputs(varargin{:});

vol = zeros(lim(1),lim(2),lim(3),'uint8');      % initialize volume
center = round([lim(1), lim(2), lim(3)] / 2);   % shape center
[X,Y,Z] = ndgrid(1:lim(1), 1:lim(2), 1:lim(3));

switch lower(varargin{1})
    case 'cube'
        vol(center - D:center + D, center - D:center + D, ...
            center - D:center + D) = 255;
    case 'cuboid'
        vol(center - D(1):center + D(1), center - D(2):center + D(2), ...
        center - D(3):center + D(3)) = 255;
    case 'sphere'
        R = sqrt((X - center(1)) .^ 2 + (Y - center(2)) .^ 2 + ...
            (Z - center(3)) .^ 2);
        vol(R <= D) = 255;
    case 'ellipsoid'    % not working yet
        R = sqrt(((X - center(1)) / D(1)) .^ 2 + ...
            ((Y - center(2)) / D(2)) .^ 2 + ...
            ((Z - center(3)) / D(3)) .^ 2);
        vol(R <= 1) = 255;
end

% Plot data
figure;
[X,Y,Z] = ind2sub(size(vol),find(vol));
scatter3(X,Y,Z,5,'*')
title(['Simulation Volume Data (',varargin{1},')'])
axis([0 lim(1) 0 lim(2) 0 lim(3)])

%%
function [D,LIM,GRADIENT] = parseinputs(varargin)
    narginchk(1,4);
    
    if ~any(strcmpi(varargin{1},{'cube','cuboid','sphere','ellipsoid'}))
        error('SHAPE should be: cube, cuboid, sphere, or ellipsoid')
    end
    
    if nargin >= 1 && nargin < 3
        switch lower(varargin{1})
            case {'cube','sphere'}
                D = 25;
            case {'cuboid','ellipsoid'}
                D = [10 20 30];
        end
        
        LIM = [100 100 100];
        GRADIENT = [0 0 0];
        
        if nargin == 2
            GRADIENT = [1 0 0];
        end
    end
    
    if nargin > 2
        if ~isnumeric(varargin{3}) || length(varargin{3}) ~= 3 || ...
            ~isequal(floor(varargin{3}), varargin{3}) || ...
            ~all(varargin{3} > 1)
            error('LIM must by a vector of length 3')
        end
        
        switch lower(varargin{1})
            case {'cube','cuboid'}
                D = varargin{2} / 2;
            case {'sphere','ellipsoid'}
                D = varargin{2};
        end
        
        LIM = varargin{3};
        GRADIENT = [0 0 0];
        
        if nargin == 4
            GRADIENT = [1 0 0];
        end
    end
    
end

end

