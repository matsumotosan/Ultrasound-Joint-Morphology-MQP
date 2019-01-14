function newPos = findPos(pose)
%UNTITLED6 Summary of this function goes here
%   function to rotate the coordinates given input angles in three axes

% Rotation angles
theta = pose(1);
phi = pose(2);
omega = pose(3);

% 3D mesh grid
dx = 5;
dy = 5;
dz = 5;
x = -50:dx:50;
y = 0:-dy:-100;
z = 0;
[X,Y,Z] = meshgrid(x,y,z);
X = X(:)'; Y = Y(:)'; Z = Z(:)';

% Plane rotations by angle around x-, y-, and z-axis
Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
Rz = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];

% Rotate coordinates
newPos = Rx * Ry * Rz * [X; Y; Z];
newPos(1,:) = newPos(1,:) + pose(4);
newPos(2,:) = newPos(2,:) + pose(5);
newPos(3,:) = newPos(3,:) + pose(6);

% figure
% scatter3(X,Y,Z);
% hold on;
% scatter3(newPos(1,:),newPos(2,:),newPos(3,:));
% xlabel('x')
% ylabel('y')
% zlabel('z')
% legend('Original','Transformed');

end

