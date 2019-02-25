function [res,noPoints] = rot_enc_res(ppr,d,rail_od,theta)
% ROT_ENC_RES Calculates spatial resolution of rotary encoder given points
% per revolution, shaft diameter, rail outer diameter, and angular range of
% view. Option to visualize resolution 
% 
% Input:       ppr - points per resolution
%                d - shaft diameter (mm)
%          rail_od - rail outer diameter (mm)
%            theta - angular range (degrees)
% 
% Output:      res - linear resolution (mm)
%         noPoints - number of points for entire path

if (theta < 0) || (theta > 360)
    error('theta must be between 0 and 360')
end

arc_len = pi * rail_od * (theta / 360); % arc length of rail
rev_len = pi * d;   % linear distance of one shaft revolution (shaft circumference)
res = rev_len / ppr;   % linear resolution
noPoints = floor(arc_len / rev_len) * ppr; % number of points in field of view

end