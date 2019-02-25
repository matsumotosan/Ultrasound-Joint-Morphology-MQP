function [intr,nonint] = findintersect(shape,volume,surface)
%FINDINTERSECT Find points on plane that intersect volume
%   
%   [INTR,NONINT] = FINDINTERSECT(SHAPE,VOLUME,SURFACE) finds the indices of
%   points in SURFACE intersecting and not intersecting with SHAPE of
%   dimensions coordinates in VOLUME. Both SURFACE and VOLUME should be
%   passed in as an n-by-3 matrix (dimensions with respect to each other
%   can differ).
%
%
%%
if size(volume,2) ~= 3
    error('VOLUME must be an n-by-3 matrix')
end

if size(surface,2) ~= 3
    error('PLANE must be an n-by-3 matrix')
end

switch lower(shape)
    case {'cube','cuboid'}
        intr = find((surface(:,1) >= min(volume(:,1)) & surface(:,1) <= max(volume(:,1))) & ...
            (surface(:,2) >= min(volume(:,2)) & surface(:,2) <= max(volume(:,2))) & ...
            (surface(:,3) >= min(volume(:,3)) & surface(:,3) <= max(volume(:,3))));
        nonint = setdiff(1:length(surface),intr);
    case {'sphere'}
        center = ceil([max(volume(:,1)) - min(volume(:,1)), ...
            max(volume(:,2)) - min(volume(:,2)), ...
            max(volume(:,3)) - min(volume(:,3))] / 2);
        r = max(rssq(vol - center), 2);
        intr = find(rssq(surface - center, 2) <= r);
        nonint = setdiff(1:length(surface),intr); 
    otherwise
        error('SHAPE must be: cube, cuboid, or sphere')
end

% Plot intersecting and non-intersecting region
% figure; hold on
% scatter3(surface(intr,1), surface(intr,2), surface(intr,3))
% scatter3(surface(nonint,1), surface(nonint,2), surface(nonint,3))
% title('Intersecting and Non-Intersecting Region of Plane')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% legend('Intersecting', 'Non-Intersecting')


end

