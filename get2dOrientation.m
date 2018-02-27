function [ output ] = get2dOrientation( p1,p2,p3 )
%GETORIENTATION Summary of this function goes here
%   Detailed explanation goes here
    val = -det([p2-p1;p3-p2]);
    output = (val>0) + -1*(val<0);
end
