function [ax,ang]=eul2aa(r)
%where r specifies rotation around the x,y,z axes, 
%ax,ang specify the same rotation as a single angle ang around an axis ax

%I'll intermediate my way through quaternions, using angle2quat:
r=deg2rad(r);
q=angle2quat(r(1),r(2),r(3),'XYZ');

if q(1)==1
    ax=[0 0 0];
    ang=0;
else
    ang = rad2deg(2 * acos(q(1)));
    ax=q(2:4);
end

