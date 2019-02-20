function mm_per_pixel = rescalc(image,depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dims = size(image);
mm_per_pixel = depth / dims;


