function pixCoord = findPos(pose,pixelCount,defaultCoord)
%findPos Find coordinates of voxels to enter pixel values given linear
%displacement and rotation data.

% Rotation angles (degrees)
theta = pose(2);    % x-axis rotation
phi = pose(3);      % y-axis rotation
omega = pose(4);    % z-axis rotation

% Frame plane (before transformation)
x = linspace(defaultCoord(1),defaultCoord(2),pixelCount(1));
y = linspace(defaultCoord(3),defaultCoord(4),pixelCount(2));
z = defaultCoord(5);
[X,Y,Z] = meshgrid(x,y,z);
X = X(:)'; Y = Y(:)'; Z = Z(:)';

% Plane rotations by angle around x-, y-, and z-axis
Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
Rz = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];

% Rotate coordinates
pixCoord = Rx * Ry * Rz * [X; Y; Z];
pixCoord(1,:) = pixCoord(1,:) + pose(5);
pixCoord(2,:) = pixCoord(2,:) + pose(6);
pixCoord(3,:) = pixCoord(3,:) + pose(7);
pixCoord = pixCoord';

figure
scatter3(X,Y,Z);
hold on;
scatter3(pixCoord(:,1),pixCoord(:,2),pixCoord(:,3));
xlabel('x')
ylabel('y')
zlabel('z')
title('Transformed Scan Plane');
legend('Original','Transformed');

end

